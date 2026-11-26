"""Assignment of reads to the best hit among the gene family.

Input is DIAMOND blastx tabular output with the columns
qseqid sseqid pident length bitscore [sstart send]. The two position columns are
optional; with them the coverage along the candidate protein is reported.

Two rules are available. Both ignore alignments shorter than ``min_aln`` aa.

symmetric (default from version 0.2.0)
  * each read goes to the gene of its best hit (highest bitscore);
  * a read of the candidate counts only if it beats the best hit to any other gene
    by at least ``margin`` bits and has at least ``min_identity`` percent identity;
  * a read of a control counts only if it beats the best hit to the candidate by at
    least ``margin`` bits; against the other genes the best hit is enough;
  * a tie never counts for the candidate or a control, whatever the gene names; a
    tie between the candidate and a control is recorded as discarded on both sides.
  Reads that cannot tell the candidate from a control are thus removed from both
  sides, so the ratio is not biased when they are nearly identical (homeologs,
  recent duplicates). The controls keep the reads they share with their own
  orthologs or paralogs labelled other, as in the census.

census (the rule of the doublecortin census, version 0.1.0)
  * each read goes to the gene of its best hit, ties broken by gene label;
  * only a read whose best hit is the candidate has to pass the margin and the
    identity; reads of the controls count without them.
  This rule reproduces the published census. It underestimates the ratio when the
  candidate and a control are nearly identical, because reads that do not tell them
  apart count for the control and not for the candidate.

The unit is the query identifier. fastq-dump with --split-spot writes both mates
of a spot under the same identifier, so the two mates of a pair are pooled and
count once, as in the original census.
"""
from collections import Counter, defaultdict
from dataclasses import dataclass, field
import hashlib

from . import diagnostics, stats

DISCARDED = "_candidate_discarded"


SYMMETRIC = "symmetric"
CENSUS = "census"


@dataclass
class Rule:
    min_aln: int = 30
    margin: float = 5.0
    min_identity: float = 75.0
    mode: str = SYMMETRIC

    def __post_init__(self):
        if self.mode not in (SYMMETRIC, CENSUS):
            raise ValueError(f"unknown rule {self.mode!r}; use {SYMMETRIC!r} or {CENSUS!r}")


@dataclass
class Counts:
    by_gene: Counter = field(default_factory=Counter)
    candidate_identity: list = field(default_factory=list)
    identity: dict = field(default_factory=lambda: defaultdict(list))     # gene -> identities of assigned reads
    candidate_spans: list = field(default_factory=list)                   # (start, end) on the candidate protein
    discarded_identity: list = field(default_factory=list)                # identity of discarded candidate reads
    control_discarded: Counter = field(default_factory=Counter)           # control -> reads discarded (symmetric)

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


def _best_per_gene(h):
    """Best hit of each gene for one read: {gene: (bitscore, pid, sstart, send)}."""
    per = {}
    for bs, g, pid, st, en in h:
        if g not in per or (bs, pid) > (per[g][0], per[g][1]):
            per[g] = (bs, pid, st, en)
    return per


def _assign_symmetric(hits, family, rule, reads):
    c = Counts()
    cand, ctrls = family.candidate, set(family.controls)
    for q, h in hits.items():
        if reads is not None and q not in reads:
            continue
        per = _best_per_gene(h)
        top_score = max(v[0] for v in per.values())
        tied = [g for g, v in per.items() if v[0] == top_score]
        if len(tied) > 1:
            tc = sorted(g for g in tied if g in ctrls)
            if cand in tied and tc:                       # candidate and control cannot be told apart
                c.by_gene[DISCARDED] += 1
                c.discarded_identity.append(per[cand][1])
                for g in tc:
                    c.control_discarded[g] += 1
                continue
            # any other tie is ambiguous and counts for neither the candidate nor a control
            others = sorted(g for g in tied if g != cand and g not in ctrls)
            if others:
                c.by_gene[others[0]] += 1
            continue
        else:
            gene = tied[0]
        bs, pid, st, en = per[gene]
        if gene == cand:
            second = max((v[0] for g, v in per.items() if g != cand), default=0.0)
            if bs - second < rule.margin or pid < rule.min_identity:
                c.by_gene[DISCARDED] += 1
                c.discarded_identity.append(pid)
                continue
            c.candidate_identity.append(pid)
            if st is not None:
                c.candidate_spans.append((min(st, en), max(st, en)))
        elif gene in ctrls:
            cand_score = per[cand][0] if cand in per else 0.0
            if bs - cand_score < rule.margin:             # cannot be told from the candidate
                c.control_discarded[gene] += 1
                continue
        c.by_gene[gene] += 1
        c.identity[gene].append(pid)
    return c


def assign(hits, family, rule=Rule(), reads=None):
    """Count reads per gene. ``reads`` restricts counting to a subset of read ids."""
    if rule.mode == SYMMETRIC:
        return _assign_symmetric(hits, family, rule, reads)
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


def summarise(counts, family, th=stats.Thresholds(), alpha=0.05, rule=None):
    ctrl = [counts.by_gene.get(g, 0) for g in family.controls]
    eff = [family.eff_length(g) for g in family.controls]
    if sum(eff) > 0:
        exp = stats.expected_reads(ctrl, eff, family.eff_length(family.candidate))
    else:
        exp = 0.0                           # no control read can be assigned: nothing to compare with
    res = stats.evaluate(counts.by_gene.get(family.candidate, 0), exp, th, alpha)
    ac = family.assignable.get(family.candidate)
    if ac is not None and ac < diagnostics.ASSIGNABLE_MIN:
        res.status = stats.UNINFORMATIVE    # the rule cannot see the candidate, so no call is possible
    out = {"candidate": family.candidate}
    if rule is not None:
        out["rule"] = rule.mode
    out.update({f"reads_{g}": n for g, n in zip(family.controls, ctrl)})
    out["candidate_discarded"] = counts.discarded
    out["control_discarded"] = sum(counts.control_discarded.values())
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
