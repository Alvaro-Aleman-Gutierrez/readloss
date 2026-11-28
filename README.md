# readloss

Test whether a gene is lost from a genome, or only missing from its assembly, using raw
whole-genome reads. Reads are assigned to their best hit among every member of the gene family,
and the count for the candidate gene is compared with the count expected from single-copy
paralogs of the same family.

    expected = (control reads / summed control length in aa) x candidate length
    ratio    = observed / expected

The ratio comes with exact Poisson limits (Garwood 1936), exact limits that also include the
uncertainty of the control count, the probability of the observed count if the gene were present,
the most likely copy-number state, and a status. By default a ratio of 0.10 or less is scored as absence,
0.30 or more as presence, values in between as intermediate, and runs with fewer than 10 expected
reads as uninformative. Intermediate runs are repeated automatically with more reads.

The method was developed and validated in a census of doublecortin across 243 specimens of
225 cypriniform species (Alemán-Gutiérrez, preprint on EcoEvoRxiv,
https://doi.org/10.32942/X2J97B).

A status of present means that reads of the coding sequence were found. It does not rule out
pseudogenization, since a gene that has lost its function can still leave reads.

## Install

readloss can be installed with conda from `environment.yml` or built as a Docker container
from the `Dockerfile`. Both are tested by continuous integration. Use the libmamba solver (the
default from conda 23.10) or mamba: the classic solver of older conda versions can run out of
memory on the bioconda dependencies.

    conda env create --solver libmamba -f environment.yml      # or: mamba env create -f environment.yml
    conda activate readloss
    pytest            # 586 tests, including the published census values

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
otherwise makes fragments of a missing paralog count for the candidate. Add the proteins of the
closest relative with an assembly, so that paralogs specific to the lineage are represented.

Two optional columns of `genes.tsv` describe the genome being tested:

- `copies`: copies per haploid genome expected in the lineage [1]. After a genome duplication
  (teleosts, salmonids, carp and goldfish, many plants) a control may be kept in two copies; say
  so, or the expected count is doubled and the candidate looks half present. When a batch mixes
  lineages with different copy numbers, give them per run in a `copies` column of `runs.tsv`
  (`dclk1a=2,dclk2a=2`) or with `--copies` for a single test.
- `location`: `autosome`, `X`, `Y`, `Z` or `W` [autosome]. Sex-linked genes are scaled by their
  dosage in the sex of the specimen, given with `--sex`, in a `sex` column of `runs.tsv`, read
  from the SRA metadata, or inferred from markers.

A gene with role `marker` and location X or Z, present in both sexes, tells the sex of a specimen
whose sex was not reported: its reads relative to the autosomal controls are about 1 in the
homogametic sex and 0.5 in the heterogametic sex.

## Use

    readloss makedb --family family.faa --out fam
    readloss test  --run SRR11971489 --db fam.dmnd --genes genes.tsv --proteins proteins.tsv
    readloss batch --runs runs.tsv   --db fam.dmnd --genes genes.tsv --proteins proteins.tsv --out results.tsv
    readloss subsample --hits run.hits.tsv --genes genes.tsv --proteins proteins.tsv
    readloss calibrate --samples known.tsv --genes genes.tsv --proteins proteins.tsv

`test` also takes `--reads` (local FASTA) or `--hits` (an existing DIAMOND table with columns
qseqid sseqid pident length bitscore). With `--stream` the download goes straight into DIAMOND
through a named pipe, so no reads are written to disk.

## What is checked before and during a run

- **Library type.** The SRA metadata of each run are read first, including whether the library
  was made without PCR (`pcr_free`: yes, no or unknown, from the library description), which makes
  GC bias much weaker. Runs that are not random
  whole-genome libraries (RNA-seq, exome or UCE capture, RAD-seq, amplicons) are scored
  uninformative with the flag `library_not_wgs`, because a gene that was not captured would look
  absent. Pools of several individuals are flagged `pooled_sample`. Use `--metadata skip` offline.
- **Aligned submissions.** Some runs (for example the 1000 Genomes high-coverage set) were
  submitted as BAM or CRAM and are stored sorted by position, so the first spots, or a few blocks,
  come from a few chromosomes only. Such runs are detected from the SRA record, are not downloaded
  unless `--spots` covers the whole run, and are scored uninformative with `aligned_sorted_run`.
- **Short reads.** When the SRA reports reads shorter than 100 nt (museum or ancient DNA), the
  minimum alignment drops to 60% of the read length in codons (at least 15 aa) unless `--min-aln`
  is given; `min_aln_used` reports the value. Rerun `selftest --length` with that read length.
- **Where the spots come from.** `--blocks 4` takes the spots in four blocks spread along the run
  instead of from its start, closer to a random sample when a run merges lanes or libraries.
- **Read sequences.** Hit tables keep the translated alignment and the whole read, so each result
  reports stop codons in candidate and control reads, the GC and the length of the reads. A
  DIAMOND that does not support these columns is detected once and the seven basic columns are
  used, leaving those diagnostics blank.

## Panel of single-copy genes and GC

A panel of unrelated single-copy genes gives a second, independent expectation and the relation
between GC and coverage in each run:

    readloss panel --faa busco_ancestral.faa --n 200 --variant-sep _ --out panel.faa
    readloss makedb --family family.faa --panel panel.faa --out fam      # writes fam.panel.tsv
    readloss test --run SRR... --db fam.dmnd --genes genes.tsv --proteins proteins.tsv \
        --panel fam.panel.tsv --cds cds.fna

- `expected_panel` and `control_vs_panel` compare the family controls with the panel. Controls
  in two copies give about 2 (flag `controls_copy_number`), and a status that changes with the
  panel expectation is flagged `panel_disagrees`.
- With `--cds` (coding sequences of candidate and controls, ids as gene or protein ids), a large
  GC difference is flagged `gc_mismatch`. Illumina libraries made with PCR under-cover GC-rich
  DNA, the main cause of genes missing from short-read assemblies, and the same bias lowers the
  reads of a GC-rich candidate. `gc_factor` is the coverage the panel predicts at the GC of the
  candidate over that at the GC of the controls; outside 0.7 to 1.4 it is flagged `gc_bias`, and
  `--gc-correct` applies it to the expected count.

The panel makes the search slower in proportion to its size; 100 to 200 genes are enough.

## Statistics

- `ci_low`, `ci_high`: exact Poisson limits of the candidate count over the expected count, as in
  the published census.
- `ci_low_cond`, `ci_high_cond`, `p_cond_if_present`: exact limits of the ratio of the two read
  rates conditional on the total, so the uncertainty of the control count is included. With few
  control reads these are noticeably wider.
- `copy_state`, `copy_state_log10_lr`: the ratio state, among 0 (absent), 0.5 (one copy where
  the controls have two, a heterozygous deletion or a hemizygous gene), 1 and 2 (duplicated),
  with the highest likelihood, and how much more likely it is than the second best, in log10 units.

## Nucleotide mode

Homeologs of an allopolyploid and recent duplicates can encode the same protein and differ at a
few nucleotides only, so a protein search cannot tell them apart. In nucleotide mode the family
is a FASTA of coding sequences (ids as in `proteins.tsv`) and reads are aligned with bwa-mem2
keeping all alignments:

    readloss makedb --family cds.fna --nt --out fam                     # bwa-mem2 index
    readloss selftest --family cds.fna --nt --genes genes.tsv --proteins proteins.tsv --out st.tsv
    readloss test --run SRR... --nt --db fam.fna --genes genes.tsv --proteins proteins.tsv \
        --assignable st.tsv
    readloss test --sam aligned.bam --genes genes.tsv --proteins proteins.tsv   # own alignments

The margin is in alignment-score points (with bwa-mem2 defaults one mismatch costs 5, so the
default margin of 5 asks for one difference; with minimap2 `-x sr` use 10), and the default
minimum identity is 95%. `selftest --nt` simulates reads of each coding sequence with sequencing
errors (`--error`) and, if asked, divergence. On simulated homeologs with both copies present,
the symmetric rule gave a ratio of 1.00 at 1.5% and at 0.5% divergence, and at 0.2%, where too
few reads differ, it returned uninformative; the census rule gave 0.14 and 0.00 (a false absence).

## Sensitive mode

    readloss sensitive --run SRR11971489 --queries queries.faa --family family.faa \
        --genes genes.tsv --proteins proteins.tsv

tBLASTn of the candidate and control proteins against the reads, then BLASTx of every read
with a hit against the whole family. It finds more divergent reads and suits distant queries
or runs with few reads. It counts reads rather than spots, applies no alignment filters, and
counts a tie between the candidate and another member for the candidate, which is
conservative for a test of absence. The genome tests of the doublecortin study used this mode with
`--rule census`, which reproduces their counts exactly; the default symmetric rule applies the
margin and identity of the DIAMOND mode, so paralog fragments no longer count for the candidate.
`--assignable` and `--copies` work here too; the fractions of `selftest` describe the DIAMOND
search, so for this mode they are an approximation.

## Worked example

`ejemplos/dcx/` holds the doublecortin family (79 proteins), the gene tables and four DIAMOND
hit tables of species of known status, so every command except the download can be tried
offline:

    readloss test --hits <(zcat ejemplos/dcx/hits/Paedocypris_micromegethes_SRR11971489.hits.tsv.gz) \
        --genes ejemplos/dcx/genes.tsv --proteins ejemplos/dcx/proteins.tsv

## Choosing the depth

The reads needed depend on genome size and read length, and in the doublecortin example they
range from 0.5 to 4.9 expected candidate reads per million spots. Estimate the yield of a run
with a short streamed pilot and download what reaches the minimum of expected reads; put the
recommended number in a `spots` column of `runs.tsv` and `batch` uses it for that run:

    readloss depth --run SRR23104325 --db fam.dmnd --genes genes.tsv --proteins proteins.tsv

With too few reads the test returns `uninformative` rather than a wrong call. The recommendation
aims at a call (10 expected reads, times the safety factor). To also exclude a partial presence,
for example to bring the upper limit of an absence below 0.10, aim at about 40 expected reads
(`readloss depth --min-expected 25`, with the default safety factor; keep the default in `test` and `batch`). `depth` also reads the SRA metadata and
gives no recommendation for runs that are not whole-genome libraries.

## Diagnostics and warnings

Every result carries diagnostics computed from the same hit table, and a `flags` column that is
empty when nothing looks wrong.

| Column | Meaning |
|---|---|
| `control_fold`, `control_p` | Agreement between controls in reads per aa, and a chi-square test |
| `ratio_loo_min`, `ratio_loo_max` | Ratio recomputed leaving out each control |
| `p_miss_one_copy` | Chance of an absence call if the candidate were present in one copy |
| `discarded_fraction`, `discarded_identity` | Candidate reads removed by the rule, and how similar they were |
| `identity_gap` | Identity of control reads that pass the candidate threshold, minus candidate identity |
| `expected_if_absent`, `p_leakage` | Candidate reads expected from control reads the rule misassigns (selftest), and the chance of the observed reads from that alone |
| `breadth`, `breadth_expected`, `largest_gap_aa` | Coverage along the candidate protein |

| Flag | What to do |
|---|---|
| `controls_disagree` | One control may not be single copy or may collect foreign reads; check it or replace it |
| `call_depends_on_one_control` | Treat the call with caution; add a control |
| `low_power_for_one_copy` | Download more spots before trusting an absence |
| `rule_may_discard_true_reads` | Close paralogs; calibrate the margin and identity for this family |
| `divergent_candidate` | Use queries from a closer species or the sensitive mode |
| `possible_partial_loss` | Part of the gene may be missing; inspect where the reads fall |
| `low_assignable_fraction` | The rule can assign under half of the candidate reads (from `selftest`); the status is set to uninformative. Use the nucleotide mode |
| `leakage_explains_candidate` | The candidate reads are no more than control reads misassigned to it (from `selftest`); the status is set to uninformative |
| `possible_pseudogene` | Many candidate reads carry stop codons: a pseudogene or a processed retrocopy may give reads although the gene is lost |
| `possible_missing_paralog` | Candidate reads beat every known member yet look like another gene, bunched in one domain; add the missing paralog |
| `gc_mismatch`, `gc_bias` | Candidate and controls differ in GC; use controls of similar GC, a panel and `--gc-correct` |
| `controls_copy_number` | Controls collect reads at a rate unlike the panel; set `copies` or choose other controls |
| `panel_disagrees` | The panel gives another status; check the controls |
| `short_reads` | Reads are short for 30 aa alignments (museum or ancient DNA); lower `--min-aln` and rerun `selftest --length` |
| `sex_unknown_for_sex_linked_gene` | Give the sex or add a marker |
| `not_expected_in_this_sex` | A Y or W gene in the homogametic sex: nothing to test |
| `library_not_wgs` | Not a random whole-genome library; the status is uninformative |
| `aligned_sorted_run` | Submitted as reads sorted by position; test the whole run or the original FASTQ |
| `pooled_sample` | Several individuals in one run; the ratio is a mixture |
| `possible_contamination` | A few candidate reads far closer to the queries than the control reads: DNA of another species (contamination, index hopping, reference DNA) |

`readloss report --results results.tsv --plot ratios.pdf` summarises calls and warnings for a batch.

## Figures to decide what to trust

`readloss plot` draws one figure per decision and writes, next to each image, a TSV with the values
drawn, so it can be redrawn in R or Python.

| Figure | Question | Command |
|---|---|---|
| `selftest` | Before a screen: can the rule assign each gene's reads to it, and where do the rest go? | `readloss plot selftest --selftest selftest.tsv --out selftest.png` |
| `power` | After a screen: which calls rest on too few expected reads? | `readloss plot power --results results.tsv --out power.png` |
| `coverage` | One run: do the candidate reads cover the whole protein, or bunch in one domain? | `readloss plot coverage --results results.tsv --run SRR... --genes genes.tsv --proteins proteins.tsv --family-fasta family.faa --out cov.png` |
| `margins` | One run: how many reads sit near the margin and identity limits of the rule? | as `coverage`, with `margins` |
| `identity` | One run: are candidate reads as similar to the queries as control reads? | as `coverage`, with `identity` |

The family includes orthologs of several species with different lengths, so read positions are
only comparable after mapping them to one reference protein per gene (the one whose length is
closest to `length_aa`). Give `--family-fasta family.faa` to `test`, `batch` and `plot`: coverage,
`breadth`, `largest_gap_aa` and `possible_partial_loss` are then exact; without it they are approximate.

A note from the example family: about a third of the control reads have below 60% identity and fall
in the kinase domain of dclk1a and dclk2a. The kinase domain is shared with many kinases that are not
in the family, so part of those reads probably come from them; part are true reads of the fast regions
of the controls in distant species. `plot margins` and `plot identity` show this mixture. Adding the
closest kinases to the family as `other` genes would absorb the foreign reads.

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

and a candidate with an assignable fraction below 0.5 is scored uninformative. The threshold
comes from real data (see Validation): X-Y gametologs at 79% protein identity (selftest 0.71)
were called correctly, while at 88-96% (0.34 and 0.37) protein mode gave false presences in
females and false absences in males. In the doublecortin family the fractions are 0.98 for DCX
and 0.90 and 0.96 for the controls. `selftest` also records the fraction of each control's
fragments given to the candidate; the result then reports `expected_if_absent`, the candidate
reads that misassigned control reads alone would produce, and a presence that this leakage
explains (`p_leakage` above 0.05) is scored uninformative. Reads simulated in another way, for example at the nucleotide level for
homeologs, can be evaluated with `selftest --hits`, naming each read `GENE::anything`.

A gene on its way to loss often evolves faster than its controls. `--divergence-table` (gene,
divergence) gives each gene its own rate, estimated for example from the identity between
orthologs in related assemblies, and `--gamma-alpha 0.5` makes rates vary along the protein in
blocks of 20 sites, as conserved domains alternate with fast linkers.

## Validation on genes of known status

Three tests with a known answer, run blind (every specimen tested under the hypothesis that the gene
is there) on public whole-genome reads, 15 to 30 million spots per run in four blocks:

- **Y gametologs of placental mammals** (truth: the sex of each specimen). DDX3Y, UTY and ZFY are on
  the Y of all placentals studied (Martínez-Pacheco et al. 2020, Genome Biol Evol 12:2015); each has
  an X gametolog 78-96% identical in protein. 26 runs of human, cattle, pig, sheep, dog and horse,
  with the X gametolog and an autosomal paralog (ZNF711 for ZFY) as controls. Horse has no Y
  sequences in RefSeq and was tested with those of the other species.
- **CHD1W of passerine birds** (truth: females ZW carry it, males ZZ do not), 8 runs of 4 species,
  CHD2 and CHD7 as controls and CHD1Z as a listed paralog.
- **PON1 of marine mammals**, inactivated by point lesions but not deleted in cetaceans, sirenians
  and two seals (Meyer et al. 2018, Science 361:591), and intact in walrus, fur seal, hippopotamus
  and five domestic species; 17 runs, PON2 and PON3 as controls. The right call is "present".

Calls (correct / wrong / uninformative), leaving out one sheep recorded as male that has no reads
of any Y gene in either mode and the X-to-autosome read ratio of a female:

| Gene | Identity to closest paralog | Protein mode | Nucleotide mode |
|---|---|---|---|
| UTY | 78-81% | 21 / 0 / 4 | 22 / 1 / 2 |
| DDX3Y | 90-94% | 0 / 0 / 25 | 21 / 1 / 3 |
| ZFY | 88-96% | 0 / 0 / 25 | 22 / 0 / 3 |
| CHD1W | close (selftest 0.34) | 0 / 0 / 8 | 8 / 0 / 0 |
| PON1 | ~60% | 11 / 0 / 6 | |

- Protein mode is right whenever it makes a call. For gametologs above ~90% identity selftest
  gives an assignable fraction below 0.5 and every call is uninformative; nucleotide mode then
  resolves them: in the 43 tests of a specimen lacking the gene, 37 had zero candidate reads and 42
  at most two.
- The two wrong calls are intermediate ratios for DDX3Y and UTY in one bull whose ZFY ratio is
  1.7; the cause is unknown.
- Nucleotide mode needs coding sequences of the same species: male horse reads matched the Y
  sequences of other species at 95.6-98.0%, lost about half their reads to the 95% threshold and
  gave false intermediate ratios, now scored uninformative (`candidate_reference_divergent`).
- PON1 was never called absent. Point lesions are rarely inside a read, so `possible_pseudogene`
  fired in one of nine inactivated genes (dugong): presence of reads says nothing about function.
- Six further runs were refused before download: four human runs stored as sorted CRAM
  (`aligned_sorted_run`) and two pig runs with a PCR library selection (`library_not_wgs`).
- Five rules came out of this validation: the assignable threshold of 0.5, the leakage
  expectation (`expected_if_absent`), the identity of control reads filtered like candidate reads,
  presence calls from an unlisted paralog (`presence_from_unlisted_paralog`, human UTY in women in
  protein mode) and the divergent nucleotide reference.

## Interpreting a call

- A call is about the specimen sequenced, not the species. Gene presence and absence vary within
  species, a heterozygous deletion gives about 0.5, and a pooled run gives a mixture. Test more
  than one individual before calling a loss for a species.
- Absent means that the coding sequence is missing from the reads. Present means that reads of
  it were found; a pseudogene or a processed retrocopy can give them, so check
  `possible_pseudogene` and, for function, use assembly-based tests of inactivating mutations.
- The test needs whole-genome reads. A gene can also be missed when it is far more divergent than
  its queries; use queries from the closest species and `selftest` with that divergence.

The counting unit is the spot. `fastq-dump --split-spot` writes both mates of a pair under
the same identifier, so a pair whose two mates hit the family counts once.

## Citation

If you use readloss, please cite the doublecortin study:
Alemán-Gutiérrez, A. (2026). Cypriniform fishes lost doublecortin repeatedly despite strong
purifying selection on the gene. EcoEvoRxiv. https://doi.org/10.32942/X2J97B

## License

MIT
