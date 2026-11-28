"""Gene family definition.

Two tab-separated files describe the family.

genes.tsv (one row per gene, with header)
    gene        gene label used in the output
    role        candidate, control, marker or other
    length_aa   protein length used for normalisation (candidate, controls and markers)
    copies      optional; copies per haploid genome expected in the lineage tested [1].
                Use 2 for a control retained in two copies after a genome duplication,
                so that the expected count is not doubled by mistake.
    location    optional; autosome, X, Y, Z or W [autosome]. Sex-linked genes are
                scaled by their dosage in the sex of the specimen (for example a Z-linked
                gene has half the dosage of an autosomal gene in a ZW female).

    A marker is a gene of known sex-linked location (X or Z) present in both sexes. Its
    reads relative to the autosomal controls tell the sex of the specimen when it is not
    given (about 1 in the homogametic sex, 0.5 in the heterogametic sex).

proteins.tsv (one row per protein in the DIAMOND database, with header)
    protein_id  identifier exactly as in the FASTA used to build the database
    gene        gene label from genes.tsv

Every protein in the database must be listed, including all paralogs that are
neither candidate nor control. Reads are assigned to the best hit among the whole
family, so a missing paralog makes its reads count for whichever member is closest.

A panel of unrelated single-copy genes (``readloss makedb --panel``) can be added to the
database; its proteins are listed in panel.tsv (protein_id, gene, length_aa) and loaded
with ``load_panel``. Panel genes compete as other genes in the assignment.
"""
import csv
from dataclasses import dataclass, field

ROLES = {"candidate", "control", "marker", "other"}
LOCATIONS = {"autosome", "X", "Y", "Z", "W"}
SEXES = {"female", "male", "unknown"}

# dosage relative to an autosomal gene in a diploid, by location and sex
DOSAGE = {
    "autosome": {"female": 1.0, "male": 1.0, "unknown": 1.0},
    "X": {"female": 1.0, "male": 0.5, "unknown": 1.0},
    "Z": {"female": 0.5, "male": 1.0, "unknown": 1.0},
    "Y": {"female": 0.0, "male": 0.5, "unknown": 0.5},
    "W": {"female": 0.5, "male": 0.0, "unknown": 0.5},
}


def norm_sex(s):
    s = (s or "").strip().lower()
    if s in ("f", "female", "woman", "hen"):
        return "female"
    if s in ("m", "male", "man", "cock"):
        return "male"
    return "unknown"


@dataclass
class Family:
    candidate: str
    controls: list
    lengths: dict                       # gene -> aa (candidate, controls and markers)
    protein_gene: dict = field(default_factory=dict)
    assignable: dict = field(default_factory=dict)   # gene -> fraction of its reads the rule can assign (selftest)
    leak: dict = field(default_factory=dict)         # control -> fraction of its fragments given to the candidate
    copies: dict = field(default_factory=dict)       # gene -> copies per haploid genome [1]
    location: dict = field(default_factory=dict)     # gene -> autosome, X, Y, Z or W
    markers: list = field(default_factory=list)      # sex-linked genes used to infer the sex
    panel: dict = field(default_factory=dict)        # panel gene -> length (aa)
    coords: dict = field(default_factory=dict)       # protein -> (reference protein, matching blocks)
    reference: dict = field(default_factory=dict)    # gene -> reference protein

    def to_reference(self, protein_id, pos):
        """Position on a family protein -> position on its gene's reference protein (1-based).

        The family holds orthologs of several species with different lengths, so a read's
        position on the protein it hit is not comparable between reads. Each protein was
        aligned with its gene's reference (``load_coords``); positions inside an aligned block
        map exactly and positions between blocks are interpolated.
        """
        m = self.coords.get(protein_id)
        if m is None or pos is None:
            return pos
        return max(1, min(m[2], self._map(m[1], pos)))

    @staticmethod
    def _map(blocks, pos):
        p = pos - 1
        prev = None
        for a, b, n in blocks:
            if a <= p < a + n:
                return b + (p - a) + 1
            if p < a:
                if prev is None:
                    return max(1, b - (a - p) + 1)
                pa, pb = prev
                frac = (p - pa) / (a - pa) if a > pa else 0
                return int(round(pb + frac * (b - pb))) + 1
            prev = (a + n, b + n)
        if prev is None:
            return pos
        return prev[1] + (p - prev[0]) + 1

    def dosage(self, gene, sex="unknown"):
        return DOSAGE[self.location.get(gene, "autosome")][norm_sex(sex)]

    def eff_length(self, gene, sex="unknown"):
        """Length that collects assignable reads: protein length x assignable fraction x copies x dosage."""
        return (self.lengths[gene] * self.assignable.get(gene, 1.0) * self.copies.get(gene, 1.0)
                * self.dosage(gene, sex))

    def sex_linked(self, gene):
        return self.location.get(gene, "autosome") != "autosome"

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
    cand, ctrls, markers, lengths, copies, loc = [], [], [], {}, {}, {}
    for r in genes:
        role = r["role"].strip().lower()
        if role not in ROLES:
            raise ValueError(f"unknown role {role!r} for gene {r['gene']!r}")
        g = r["gene"]
        if role == "candidate":
            cand.append(g)
        elif role == "control":
            ctrls.append(g)
        elif role == "marker":
            markers.append(g)
        if role in ("candidate", "control", "marker"):
            lengths[g] = float(r["length_aa"])
        if (r.get("copies") or "").strip():
            copies[g] = float(r["copies"])
            if copies[g] <= 0:
                raise ValueError(f"copies of {g!r} must be positive")
        l_ = (r.get("location") or "").strip() or "autosome"
        l_ = l_ if l_ == "autosome" else l_.upper()
        if l_ not in LOCATIONS:
            raise ValueError(f"unknown location {l_!r} for gene {g!r}; use autosome, X, Y, Z or W")
        if l_ != "autosome":
            loc[g] = l_
    if len(cand) != 1:
        raise ValueError(f"genes.tsv needs exactly one candidate, found {len(cand)}")
    if not ctrls:
        raise ValueError("genes.tsv needs at least one control")
    for m in markers:
        if loc.get(m) not in ("X", "Z"):
            raise ValueError(f"marker {m!r} needs location X or Z")
    known = {r["gene"] for r in genes}
    pg = {}
    for r in _rows(proteins_tsv):
        if r["gene"] not in known:
            raise ValueError(f"protein {r['protein_id']!r} maps to unknown gene {r['gene']!r}")
        pg[r["protein_id"]] = r["gene"]
    return Family(cand[0], ctrls, lengths, pg, copies=copies, location=loc, markers=markers)


def load_assignable(path, family):
    """Read the assignable fractions written by ``readloss selftest`` into ``family``."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    a = {r["gene"]: float(r["assignable_fraction"]) for r in rows}
    missing = [g for g in [family.candidate, *family.controls] if g not in a]
    if missing:
        raise ValueError(f"selftest table lacks {', '.join(missing)}; run selftest with the same genes.tsv")
    family.assignable = {g: a[g] for g in [family.candidate, *family.controls]}
    # fraction of each control's fragments that the rule gave to the candidate (false candidate reads)
    family.leak = {}
    for r in rows:
        if r["gene"] in family.controls and r.get("assigned_candidate") not in (None, "") and \
                float(r.get("fragments") or 0) > 0:
            family.leak[r["gene"]] = float(r["assigned_candidate"]) / float(r["fragments"])
    return family


def load_panel(path, family):
    """Add the panel genes of ``readloss makedb --panel`` (panel.tsv) to ``family``."""
    for r in _rows(path):
        g = r["gene"]
        if g in family.lengths or g in set(family.protein_gene.values()):
            raise ValueError(f"panel gene {g!r} clashes with a family gene")
        family.protein_gene[r["protein_id"]] = g
        family.panel[g] = float(r["length_aa"])
    return family


def load_coords(family_faa, family):
    """Align every candidate, control and marker protein with its gene's reference protein.

    The reference of a gene is its protein whose length is closest to length_aa in genes.tsv,
    so that coverage and breadth are measured on one protein of that length.
    """
    import difflib
    from .selftest import read_fasta
    seqs = read_fasta(family_faa)
    genes = [family.candidate, *family.controls, *family.markers]
    by = {}
    for pid, g in family.protein_gene.items():
        if g in genes and pid in seqs:
            by.setdefault(g, []).append(pid)
    for g, pids in by.items():
        ref = min(pids, key=lambda x: (abs(len(seqs[x]) - family.lengths[g]), x))
        family.reference[g] = ref
        for pid in pids:
            sm = difflib.SequenceMatcher(None, seqs[pid], seqs[ref], autojunk=False)
            family.coords[pid] = (ref, [(bl.a, bl.b, bl.size) for bl in sm.get_matching_blocks() if bl.size],
                                  len(seqs[ref]))
    return family
