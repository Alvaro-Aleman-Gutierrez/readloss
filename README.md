# readloss

Test whether a gene is lost from a genome, or only missing from its assembly, using raw
whole-genome reads. Reads are assigned to their best hit among every member of the gene family,
and the count for the candidate gene is compared with the count expected from single-copy
paralogs of the same family.

    expected = (control reads / summed control length in aa) x candidate length
    ratio    = observed / expected

The ratio comes with exact Poisson limits (Garwood 1936), the probability of the observed count
if the gene were present, and a status. By default a ratio of 0.10 or less is scored as absence,
0.30 or more as presence, values in between as intermediate, and runs with fewer than 10 expected
reads as uninformative. Intermediate runs are repeated automatically with more reads.

The method was developed and validated in a census of doublecortin across 243 specimens of
225 cypriniform species (Alemán-Gutiérrez, preprint on EcoEvoRxiv,
https://doi.org/10.32942/X2J97B).

A status of present means that reads of the coding sequence were found. It does not rule out
pseudogenization, since a gene that has lost its function can still leave reads.

## Install

readloss can be installed with conda from `environment.yml` or built as a Docker container
from the `Dockerfile`. Both are tested by continuous integration.

    conda env create -f environment.yml
    conda activate readloss
    pytest            # 534 tests, including the published census values

or with the container

    docker build -t readloss .
    docker run --rm -v "$PWD":/data readloss --help

## Describe the family

`genes.tsv`

    gene    role       length_aa
    DCX     candidate  357
    dclk1a  control    770
    dclk2a  control    826
    dclk3   other      0

`proteins.tsv` lists every protein in the family FASTA and the gene it belongs to.
Include all paralogs, also those that are neither candidate nor control. A conserved domain
otherwise makes fragments of a missing paralog count for the candidate.

## Use

    readloss makedb --family family.faa --out fam
    readloss test  --run SRR11971489 --db fam.dmnd --genes genes.tsv --proteins proteins.tsv
    readloss batch --runs runs.tsv   --db fam.dmnd --genes genes.tsv --proteins proteins.tsv --out results.tsv
    readloss subsample --hits run.hits.tsv --genes genes.tsv --proteins proteins.tsv
    readloss calibrate --samples known.tsv --genes genes.tsv --proteins proteins.tsv

`test` also takes `--reads` (local FASTA) or `--hits` (an existing DIAMOND table with columns
qseqid sseqid pident length bitscore). With `--stream` the download goes straight into DIAMOND
through a named pipe, so no reads are written to disk.

## Sensitive mode

    readloss sensitive --run SRR11971489 --queries queries.faa --family family.faa \
        --genes genes.tsv --proteins proteins.tsv

tBLASTn of the candidate and control proteins against the reads, then BLASTx of every read
with a hit against the whole family. It finds more divergent reads and suits distant queries
or runs with few reads. It counts reads rather than spots, applies no alignment filters, and
counts a tie between the candidate and another member for the candidate, which is
conservative for a test of absence. This is the version behind the genome tests of the
doublecortin study, whose counts it reproduces exactly.

## Worked example

`ejemplos/dcx/` holds the doublecortin family (79 proteins), the gene tables and four DIAMOND
hit tables of species of known status, so every command except the download can be tried
offline:

    readloss test --hits <(zcat ejemplos/dcx/hits/Paedocypris_micromegethes_SRR11971489.hits.tsv.gz) \
        --genes ejemplos/dcx/genes.tsv --proteins ejemplos/dcx/proteins.tsv

## Choosing the depth

The reads needed depend on genome size and read length, and in the doublecortin example they
range from 0.5 to 4.9 expected candidate reads per million spots. Estimate the yield of a run
with a short streamed pilot and download what reaches the minimum of expected reads:

    readloss depth --run SRR23104325 --db fam.dmnd --genes genes.tsv --proteins proteins.tsv

With too few reads the test returns `uninformative` rather than a wrong call.

## Diagnostics and warnings

Every result carries diagnostics computed from the same hit table, and a `flags` column that is
empty when nothing looks wrong.

| Column | Meaning |
|---|---|
| `control_fold`, `control_p` | Agreement between controls in reads per aa, and a chi-square test |
| `ratio_loo_min`, `ratio_loo_max` | Ratio recomputed leaving out each control |
| `p_miss_one_copy` | Chance of an absence call if the candidate were present in one copy |
| `discarded_fraction`, `discarded_identity` | Candidate reads removed by the rule, and how similar they were |
| `identity_gap` | Control identity minus candidate identity |
| `breadth`, `breadth_expected`, `largest_gap_aa` | Coverage along the candidate protein |

| Flag | What to do |
|---|---|
| `controls_disagree` | One control may not be single copy or may collect foreign reads; check it or replace it |
| `call_depends_on_one_control` | Treat the call with caution; add a control |
| `low_power_for_one_copy` | Download more spots before trusting an absence |
| `rule_may_discard_true_reads` | Close paralogs; calibrate the margin and identity for this family |
| `divergent_candidate` | Use queries from a closer species or the sensitive mode |
| `possible_partial_loss` | Part of the gene may be missing; inspect where the reads fall |
| `low_assignable_fraction` | The rule can assign few candidate reads (from `selftest`); the status is set to uninformative |

`readloss report --results results.tsv --plot ratios.pdf` summarises calls and warnings for a batch.

## Large screens with Snakemake

    snakemake -s workflow/Snakefile --configfile workflow/config.yaml --cores 16

Each run is an independent job, so the same workflow runs on a cluster with the executor of the
scheduler (for example `--executor slurm --jobs 50`). Failed downloads are retried, and the
results are gathered into one table with the report and a figure.

## Assignment rule

Alignments of at least 30 aa. Each read goes to the gene of its best hit.

- A read of the candidate counts only if it beats the best hit to any other gene by at least
  5 bits and has at least 75% identity.
- A read of a control counts only if it beats the best hit to the candidate by at least 5 bits.
- A tie never counts for the candidate or a control, whatever the gene names.

This is the `symmetric` rule, the default from version 0.2.0. Reads that cannot tell the
candidate from a control are removed from both sides, so the ratio is not biased when the two
are nearly identical, as homeologs of an allopolyploid or recent duplicates. The rule of
version 0.1.0, under which only candidate reads had to pass the margin and ties went to the
gene with the later name, is kept as `--rule census` and reproduces the doublecortin census
exactly. On the four example species the two rules give the same status and ratios within
0.03, but with a candidate nearly identical to its control the census rule can call a present
gene absent.

The rule is set with `--rule`, `--min-aln`, `--margin` and `--min-identity` and can be
calibrated with `readloss calibrate` on species of known status.

## Detectability: selftest

The ratio assumes that reads of the candidate and of the controls can be assigned to their
genes. Run `selftest` once per family to check it:

    readloss selftest --family family.faa --genes genes.tsv --proteins proteins.tsv --out selftest.tsv
    readloss test --run SRR11971489 --db fam.dmnd --genes genes.tsv --proteins proteins.tsv \
        --assignable selftest.tsv

`selftest` cuts every candidate and control protein into fragments the length of a read
(`--length 50` aa for 150 nt reads), can add random substitutions to mimic a divergent species
(`--divergence 0.1`), searches them against the family with the settings of the real search and
applies the same rule. It reports, for each gene, the fraction of fragments assigned back to it.
With `--assignable` the expected count uses those fractions,

    expected = control reads / sum(control length x assignable) x candidate length x assignable

and a candidate with an assignable fraction below 0.2 is scored uninformative, because an absence
call would mean nothing. In the doublecortin family the fractions are 0.98 for DCX and 0.90 and
0.96 for the controls. Reads simulated in another way, for example at the nucleotide level for
homeologs, can be evaluated with `selftest --hits`, naming each read `GENE::anything`.

The counting unit is the spot. `fastq-dump --split-spot` writes both mates of a pair under
the same identifier, so a pair whose two mates hit the family counts once.

## Citation

If you use readloss, please cite the doublecortin study:
Alemán-Gutiérrez, A. (2026). Cypriniform fishes lost doublecortin repeatedly despite strong
purifying selection on the gene. EcoEvoRxiv. https://doi.org/10.32942/X2J97B

## License

MIT
