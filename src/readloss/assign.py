"""Assignment of reads to the best hit among the gene family.

Input is DIAMOND blastx tabular output with the columns
qseqid sseqid pident length bitscore [sstart send]. The two position columns are
optional; with them the coverage along the candidate protein is reported.

Rule (calibrated in the doublecortin study):
  * alignments shorter than ``min_aln`` aa are ignored;
  * each read goes to the gene of its best hit (highest bitscore);
  * a read whose best hit is the candidate counts for the candidate only if the
    best hit beats the best non-candidate hit by at least ``margin`` bits and has
    at least ``min_identity`` percent identity; otherwise it is discarded.

The unit is the query identifier. fastq-dump with --split-spot writes both mates
of a spot under the same identifier, so the two mates of a pair are pooled and
count once, as in the original census.
"""
from collections import Counter, defaultdict
from dataclasses import dataclass, field
import hashlib

from . import diagnostics, stats

DISCARDED = "_candidate_discarded"


@dataclass
class Rule:
    min_aln: int = 30
    margin: float = 5.0
    min_identity: float = 75.0


@dataclass
class Counts:
    by_gene: Counter = field(default_factory=Counter)
    candidate_identity: list = field(default_factory=list)
    identity: dict = field(default_factory=lambda: defaultdict(list))     # gene -> identities of assigned reads
    candidate_spans: list = field(default_factory=list)                   # (start, end) on the candidate protein
    discarded_identity: list = field(default_factory=list)                # identity of discarded candidate reads

    @property
    def discarded(self):
        return self.by_gene.get(DISCARDED, 0)


def read_hits(path, family, rule=Rule()):
    """Group hits by read: {read: [(bitscore, gene, pident, sstart, send), ...]}.

    Positions are None for tables with five columns. The tuple order keeps the tie
    rule of the original census: bitscore, then gene label, then identity.
    """
    hits = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.split()
            q, s, pid, ln, bs = f[:5]
            st, en = (int(f[5]), int(f[6])) if len(f) >= 7 else (None, None)
            if int(float(ln)) >= rule.min_aln:
                hits[q].append((float(bs), family.gene_of(s), float(pid), st, en))
    return hits


def assign(hits, family, rule=Rule(), reads=None):
    """Count reads per gene. ``reads`` restricts counting to a subset of read ids."""
    c = Counts()
    cand = family.candidate
    for q, h in hits.items():
        if reads is not None and q not in reads:
            continue
        best = max(h)                       # same order as sorting (bitscore, gene, pid) descending
        if best[1] == cand:
            other = max((x[0] for x in h if x[1] != cand), default=0)
            if best[0] - other < rule.margin or best[2] < rule.min_identity:
                c.by_gene[DISCARDED] += 1
                c.discarded_identity.append(best[2])
                continue
            c.candidate_identity.append(best[2])
            if best[3] is not None:
                c.candidate_spans.append((min(best[3], best[4]), max(best[3], best[4])))
        c.by_gene[best[1]] += 1
        c.identity[best[1]].append(best[2])
    return c


def summarise(counts, family, th=stats.Thresholds(), alpha=0.05):
    ctrl = [counts.by_gene.get(g, 0) for g in family.controls]
    exp = stats.expected_reads(ctrl, [family.lengths[g] for g in family.controls],
                               family.lengths[family.candidate])
    res = stats.evaluate(counts.by_gene.get(family.candidate, 0), exp, th, alpha)
    out = {"candidate": family.candidate}
    out.update({f"reads_{g}": n for g, n in zip(family.controls, ctrl)})
    out["candidate_discarded"] = counts.discarded
    out.update(res.as_dict())
    ids = counts.candidate_identity
    out["mean_identity_candidate"] = sum(ids) / len(ids) if ids else 0.0
    out.update(diagnostics.run(counts, family, res, th))
    return out


def subsample_ids(hits, fraction, seed=1):
    """Deterministic thinning of reads: keeps each read with probability ``fraction``.

    Reads without hits never enter the table, but thinning all spots uniformly
    thins the reads with hits in the same way, so a hit table can be subsampled
    without downloading the reads again.
    """
    lim = int(fraction * 2**64)
    salt = str(seed).encode()
    return {q for q in hits
            if int.from_bytes(hashlib.blake2b(q.encode(), digest_size=8, key=salt).digest(), "big") < lim}
