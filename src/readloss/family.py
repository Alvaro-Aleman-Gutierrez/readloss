"""Gene family definition.

Two tab-separated files describe the family.

genes.tsv (one row per gene, with header)
    gene        gene label used in the output
    role        candidate, control or other
    length_aa   protein length used for normalisation (candidate and controls)

proteins.tsv (one row per protein in the DIAMOND database, with header)
    protein_id  identifier exactly as in the FASTA used to build the database
    gene        gene label from genes.tsv

Every protein in the database must be listed, including all paralogs that are
neither candidate nor control. Reads are assigned to the best hit among the whole
family, so a missing paralog makes its reads count for whichever member is closest.
"""
import csv
from dataclasses import dataclass, field

ROLES = {"candidate", "control", "other"}


@dataclass
class Family:
    candidate: str
    controls: list
    lengths: dict                       # gene -> aa (candidate and controls)
    protein_gene: dict = field(default_factory=dict)
    assignable: dict = field(default_factory=dict)   # gene -> fraction of its reads the rule can assign (selftest)

    def eff_length(self, gene):
        """Length that collects assignable reads: protein length times the assignable fraction (1 if unknown)."""
        return self.lengths[gene] * self.assignable.get(gene, 1.0)

    def gene_of(self, protein_id):
        try:
            return self.protein_gene[protein_id]
        except KeyError:
            raise KeyError(f"protein {protein_id!r} is not listed in proteins.tsv") from None


def _rows(path):
    with open(path, newline="") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")]
    return rows


def load(genes_tsv, proteins_tsv):
    genes = _rows(genes_tsv)
    cand, ctrls, lengths = [], [], {}
    for r in genes:
        role = r["role"].strip().lower()
        if role not in ROLES:
            raise ValueError(f"unknown role {role!r} for gene {r['gene']!r}")
        if role == "candidate":
            cand.append(r["gene"])
        elif role == "control":
            ctrls.append(r["gene"])
        if role in ("candidate", "control"):
            lengths[r["gene"]] = float(r["length_aa"])
    if len(cand) != 1:
        raise ValueError(f"genes.tsv needs exactly one candidate, found {len(cand)}")
    if not ctrls:
        raise ValueError("genes.tsv needs at least one control")
    known = {r["gene"] for r in genes}
    pg = {}
    for r in _rows(proteins_tsv):
        if r["gene"] not in known:
            raise ValueError(f"protein {r['protein_id']!r} maps to unknown gene {r['gene']!r}")
        pg[r["protein_id"]] = r["gene"]
    return Family(cand[0], ctrls, lengths, pg)


def load_assignable(path, family):
    """Read the assignable fractions written by ``readloss selftest`` into ``family``."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    a = {r["gene"]: float(r["assignable_fraction"]) for r in rows}
    missing = [g for g in [family.candidate, *family.controls] if g not in a]
    if missing:
        raise ValueError(f"selftest table lacks {', '.join(missing)}; run selftest with the same genes.tsv")
    family.assignable = {g: a[g] for g in [family.candidate, *family.controls]}
    return family
