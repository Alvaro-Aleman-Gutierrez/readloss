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

from scipy.stats import poisson

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
    discarded_spans: list = field(default_factory=list)                   # spans of discarded candidate reads
    identity_only_spans: list = field(default_factory=list)               # discarded only for low identity: the
    #                                     read beats every known member, so its source may be missing from the family
    stops: Counter = field(default_factory=Counter)                       # gene -> counted reads with a stop codon
    with_seq: Counter = field(default_factory=Counter)                    # gene -> counted reads with a sequence
    gc: dict = field(default_factory=lambda: defaultdict(list))           # gene -> GC of counted reads
    read_len: list = field(default_factory=list)                          # length (nt) of counted reads
    spans: dict = field(default_factory=lambda: defaultdict(list))        # gene -> spans of counted reads
    margins: list = field(default_factory=list)                           # (gene, margin, identity, kept) of reads
    #                                     whose best hit is the candidate or a control (margin over the competitor)

    @property
    def discarded(self):
        return self.by_gene.get(DISCARDED, 0)


def _gc(seq):
    s = seq.upper()
    acgt = sum(s.count(b) for b in "ACGT")
    return (s.count("G") + s.count("C")) / acgt if acgt else float("nan")


def read_hits(path, family, rule=Rule(), info=None):
    """Group hits by read: {read: [(bitscore, gene, pident, sstart, send), ...]}.

    Positions are None for tables with five columns. The tuple order keeps the tie
    rule of the original census: bitscore, then gene label, then identity.

    Tables written by readloss 0.2.0 or later add two columns, the translated aligned
    read (qseq_translated) and the whole read (full_qseq). If ``info`` is a dict, it is
    filled with {(read, gene): (has_stop, gc, length_nt)} for the best hit of each gene,
    used for the pseudogene, GC and read-length diagnostics.
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
                g, b = family.gene_of(s), float(bs)
                if st is not None and family.coords:
                    st, en = family.to_reference(s, st), family.to_reference(s, en)
                hits[q].append((b, g, float(pid), st, en))
                if info is not None and len(f) >= 9 and f[8] != "*":
                    k = (q, g)
                    if k not in info or b > info[k][0]:
                        info[k] = (b, "*" in f[7], _gc(f[8]), len(f[8]))
    if info is not None:
        for k, v in info.items():
            if len(v) == 4:
                info[k] = v[1:]
    return hits


def _record(c, info, q, gene):
    if not info:
        return
    x = info.get((q, gene))
    if x is None:
        return
    stop, gc, ln = x
    c.with_seq[gene] += 1
    if stop:
        c.stops[gene] += 1
    if gc == gc:
        c.gc[gene].append(gc)
    c.read_len.append(ln)


def _best_per_gene(h):
    """Best hit of each gene for one read: {gene: (bitscore, pid, sstart, send)}."""
    per = {}
    for bs, g, pid, st, en in h:
        if g not in per or (bs, pid) > (per[g][0], per[g][1]):
            per[g] = (bs, pid, st, en)
    return per


def _assign_symmetric(hits, family, rule, reads, info=None):
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
                c.margins.append((cand, 0.0, per[cand][1], False))
                for g in tc:
                    c.margins.append((g, 0.0, per[g][1], False))
                if per[cand][2] is not None:
                    c.discarded_spans.append((min(per[cand][2], per[cand][3]), max(per[cand][2], per[cand][3])))
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
            c.margins.append((cand, bs - second, pid, bs - second >= rule.margin and pid >= rule.min_identity))
            if bs - second < rule.margin or pid < rule.min_identity:
                c.by_gene[DISCARDED] += 1
                c.discarded_identity.append(pid)
                if st is not None:
                    c.discarded_spans.append((min(st, en), max(st, en)))
                    if bs - second >= rule.margin:
                        c.identity_only_spans.append((min(st, en), max(st, en)))
                continue
            c.candidate_identity.append(pid)
            if st is not None:
                c.candidate_spans.append((min(st, en), max(st, en)))
        elif gene in ctrls:
            cand_score = per[cand][0] if cand in per else 0.0
            c.margins.append((gene, bs - cand_score, pid, bs - cand_score >= rule.margin))
            if bs - cand_score < rule.margin:             # cannot be told from the candidate
                c.control_discarded[gene] += 1
                continue
        c.by_gene[gene] += 1
        c.identity[gene].append(pid)
        if st is not None:
            c.spans[gene].append((min(st, en), max(st, en)))
        _record(c, info, q, gene)
    return c


def assign(hits, family, rule=Rule(), reads=None, info=None):
    """Count reads per gene. ``reads`` restricts counting to a subset of read ids.

    ``info`` is the dict filled by read_hits; with it the counts also hold the GC, stop
    codons and length of the counted reads.
    """
    if rule.mode == SYMMETRIC:
        return _assign_symmetric(hits, family, rule, reads, info)
    c = Counts()
    cand = family.candidate
    for q, h in hits.items():
        if reads is not None and q not in reads:
            continue
        best = max(h)                       # same order as sorting (bitscore, gene, pid) descending
        if best[1] in family.controls:
            cs = max((x[0] for x in h if x[1] == cand), default=0.0)
            c.margins.append((best[1], best[0] - cs, best[2], True))
        if best[1] == cand:
            other = max((x[0] for x in h if x[1] != cand), default=0)
            c.margins.append((cand, best[0] - other, best[2], best[0] - other >= rule.margin and
                              best[2] >= rule.min_identity))
            if best[0] - other < rule.margin or best[2] < rule.min_identity:
                c.by_gene[DISCARDED] += 1
                c.discarded_identity.append(best[2])
                if best[3] is not None:
                    c.discarded_spans.append((min(best[3], best[4]), max(best[3], best[4])))
                    if best[0] - other >= rule.margin:
                        c.identity_only_spans.append((min(best[3], best[4]), max(best[3], best[4])))
                continue
            c.candidate_identity.append(best[2])
            if best[3] is not None:
                c.candidate_spans.append((min(best[3], best[4]), max(best[3], best[4])))
        c.by_gene[best[1]] += 1
        c.identity[best[1]].append(best[2])
        if best[3] is not None:
            c.spans[best[1]].append((min(best[3], best[4]), max(best[3], best[4])))
        _record(c, info, q, best[1])
    return c


def summarise(counts, family, th=stats.Thresholds(), alpha=0.05, rule=None, sex="unknown", cds_gc=None,
              gc_correct=False):
    """One result row: counts, ratio, limits, copy-number state and diagnostics.

    sex          sex of the specimen (female, male or unknown); unknown is inferred from markers
    cds_gc       {gene: GC of its coding sequence}, used for the GC diagnostics
    gc_correct   multiply the expected count by the GC factor estimated from the panel
    """
    from . import extras, family as fam
    cand = family.candidate
    sex, sex_source = fam.norm_sex(sex), "given"
    ctrl = [counts.by_gene.get(g, 0) for g in family.controls]
    marker_ratio = float("nan")
    if sex == "unknown":
        sex_source = "unknown"
        auto = [g for g in family.controls if not family.sex_linked(g)]
        al = sum(family.eff_length(g) for g in auto)
        rate = sum(counts.by_gene.get(g, 0) for g in auto) / al if al else float("nan")
        s_, marker_ratio = extras.infer_sex(counts, family, rate)
        if s_ != "unknown":
            sex, sex_source = s_, "markers"
    eff = [family.eff_length(g, sex) for g in family.controls]
    lc = family.eff_length(cand, sex)
    exp = stats.expected_reads(ctrl, eff, lc) if sum(eff) > 0 else 0.0
    # GC and panel
    gc_c = (cds_gc or {}).get(cand, extras._mean(counts.gc.get(cand, [])))
    gk = [(cds_gc or {}).get(g) for g in family.controls]
    gc_k = extras._mean([x for x in gk if x is not None]) if any(x is not None for x in gk) else \
        extras._mean([x for g in family.controls for x in counts.gc.get(g, [])])
    prate, npanel = extras.panel_rate(counts, family) if family.panel else (float("nan"), 0)
    curve = extras.gc_curve(counts, family) if family.panel else []
    fc, fk = extras.curve_at(curve, gc_c), extras.curve_at(curve, gc_k)
    gc_factor = fc / fk if fc == fc and fk == fk and fk > 0 else float("nan")
    if gc_correct and gc_factor == gc_factor:
        exp *= gc_factor
    obs = counts.by_gene.get(cand, 0)
    res = stats.evaluate(obs, exp, th, alpha)
    ac = family.assignable.get(cand)
    if ac is not None and ac < diagnostics.ASSIGNABLE_MIN:
        res.status = stats.UNINFORMATIVE    # the rule cannot see the candidate, so no call is possible
    if lc <= 0:
        res.status = stats.UNINFORMATIVE    # the candidate is not expected in this sex
    # false candidate reads expected from the controls (selftest): a presence call needs more than that
    bg, p_leak = float("nan"), float("nan")
    if family.leak:
        bg = sum(counts.by_gene.get(g, 0) / family.assignable[g] * family.leak.get(g, 0.0)
                 for g in family.controls if family.assignable.get(g, 0) > 0)
        p_leak = float(poisson.sf(obs - 1, bg)) if obs > 0 and bg > 0 else (1.0 if obs == 0 else 0.0)
    leaked = obs > 0 and p_leak == p_leak and p_leak > diagnostics.LEAK_P and \
        res.status in (stats.PRESENT, stats.INTERMEDIATE)
    if leaked:
        res.status = stats.UNINFORMATIVE    # the candidate reads could all be misassigned control reads
    lk = sum(eff) * (1 / gc_factor if gc_correct and gc_factor == gc_factor and gc_factor > 0 else 1)
    clo, chi, cp = stats.conditional(obs, sum(ctrl), lc, lk, alpha)
    cn, cn_lr, _ = stats.copy_number_states(obs, sum(ctrl), lc, lk)
    out = {"candidate": cand}
    if rule is not None:
        out["rule"] = rule.mode
    out.update({f"reads_{g}": n for g, n in zip(family.controls, ctrl)})
    out["candidate_discarded"] = counts.discarded
    out["control_discarded"] = sum(counts.control_discarded.values())
    out.update(res.as_dict())
    out.update({"ci_low_cond": clo, "ci_high_cond": chi, "p_cond_if_present": cp,
                "copy_state": cn, "copy_state_log10_lr": cn_lr})
    ids = counts.candidate_identity
    out["mean_identity_candidate"] = sum(ids) / len(ids) if ids else 0.0
    ctx = {"sex": sex, "sex_source": sex_source, "marker_ratio": marker_ratio, "gc_candidate": gc_c,
           "gc_controls": gc_k, "gc_factor": gc_factor, "panel_rate": prate, "panel_genes": npanel,
           "expected": exp, "lc": lc, "min_identity": rule.min_identity if rule is not None else 0.0}
    out.update(diagnostics.run(counts, family, res, th, ctx))
    out["expected_if_absent"] = bg if bg == bg else ""
    out["p_leakage"] = p_leak if p_leak == p_leak else ""
    fl = [x for x in out["flags"].split(",") if x]
    if leaked:
        fl.append("leakage_explains_candidate")
    # a presence carried by reads that look like an unlisted family member (bunched in one region and less
    # similar than the controls) is not evidence for the candidate: human UTY in women, protein mode
    if out["status"] == stats.PRESENT and "possible_missing_paralog" in fl and "divergent_candidate" in fl:
        out["status"] = stats.UNINFORMATIVE
        fl.append("presence_from_unlisted_paralog")
    # nucleotide mode with a candidate sequence of another species: true reads are lost, so a low ratio is not
    # evidence of absence; a presence still is
    if "candidate_reference_divergent" in fl and out["status"] in (stats.ABSENT, stats.INTERMEDIATE):
        out["status"] = stats.UNINFORMATIVE
    out["flags"] = ",".join(fl)
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
