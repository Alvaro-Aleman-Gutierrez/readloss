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
