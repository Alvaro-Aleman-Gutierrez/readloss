# Versions and dates behind the doublecortin example

## Software

| Program | Census (7 to 8 Nov 2026) | Example tables regenerated (16 Nov 2026) | Tested environment |
|---|---|---|---|
| DIAMOND | 2.1.10 | 2.1.10 | 2.1.10 |
| SRA Toolkit (fastq-dump) | 2.11.3 | 2.11.3 | 3.4.1 (same output format) |
| BLAST+ (sensitive mode) | 2.12.0 | 2.12.0 | 2.17.0 (solved, tests pass) |
| Python / SciPy | 3.11 / 1.13 | 3.11.5 / 1.13.1 | 3.11.18 / conda-forge |
| readloss | scripts it replaces | 0.1.0 | 0.1.0 |

## Family database

- `family.faa`, 79 proteins of the doublecortin family, MD5 `b9832e0a0bc003eab5a234794ea5fc08`.
- The DIAMOND database used in the census contains exactly these 79 sequences, which was checked
  with `diamond getseq`. Rebuild it with `readloss makedb --family family.faa --out fam`.

## Search settings

DIAMOND blastx, E < 1e-5, five targets per read, one HSP, block size 0.5, one index chunk,
output columns qseqid sseqid pident length bitscore. Assignment with alignments of at least
30 aa, a margin of 5 bits and at least 75% identity.

## Example hit tables

| Species | Run | Spots | Origin |
|---|---|---|---|
| *Danionella cerebrum* | ERR3284980 | 15,000,000 | census, 7 Nov 2026 |
| *Paedocypris micromegethes* | SRR11971489 | 15,000,000 | regenerated with `readloss test --stream`, 16 Nov 2026 |
| *Ictalurus furcatus* | SRR18966194 | 40,000,000 | independent control, 7 Nov 2026 |
| *Paramisgurnus dabryanus* | SRR23104325 | 35,000,000 | independent control, 7 Nov 2026 |

The calibration value of the doublecortin study for *P. micromegethes* (ratio 0.016) comes from
the 12,000,000-spot reads of its genome test. Searching those reads again on 16 Nov 2026 gave a
hit table byte-identical to the original calibration table of 7 Nov 2026, and readloss
reproduces the value exactly (1 read against 63.1 expected). With the 15,000,000 spots of the
table above the ratio is 0.014.

Without the assignment rule the same 12,000,000-spot table gives 8 candidate reads (ratio 0.127,
mean identity 68%), which are fragments of paralogs. The margin and identity filters discard 7
of them. This is the example of a bad result used in the protocol.
