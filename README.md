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

The method was developed and validated in a census of doublecortin across 243 cypriniform
specimens (Alemán-Gutiérrez, preprint in preparation).

## Install

    conda env create -f environment.yml
    conda activate readloss
    pytest            # 512 tests, including the published census values

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

## Assignment rule

Alignments of at least 30 aa. A read counts for the candidate only if its best hit is the
candidate, beats the best non-candidate hit by at least 5 bits, and has at least 75% identity.
The rule is set with `--min-aln`, `--margin` and `--min-identity` and can be calibrated with
`readloss calibrate` on species of known status.

The counting unit is the spot. `fastq-dump --split-spot` writes both mates of a pair under
the same identifier, so a pair whose two mates hit the family counts once.

## License

MIT
