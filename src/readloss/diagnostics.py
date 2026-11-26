"""Diagnostics added to every result, so that users can judge a call before trusting it.

control_fold         highest over lowest reads per aa among the controls (1 = perfect agreement)
control_p            chi-square test that the controls collect reads in proportion to their length
ratio_loo_min/max    ratio recomputed leaving out each control in turn
p_miss_one_copy      probability of a false absence call if the candidate were present in one copy
                     while the controls have two (ratio 0.5), given the expected reads
discarded_fraction   candidate reads discarded by the margin or identity rule, over all candidate reads
discarded_identity   mean identity of the discarded reads (low = paralog fragments, high = possibly true reads)
identity_gap         mean identity of control reads minus mean identity of candidate reads
breadth              fraction of the candidate protein covered by kept reads (needs sstart and send)
breadth_expected     breadth expected for that number of reads if they fell at random along the protein
largest_gap_aa       longest stretch of the candidate protein without reads
assignable_candidate fraction of candidate reads the rule can assign, from selftest (blank without it)
flags                comma-separated warnings (empty when nothing looks wrong)
"""
import math

from scipy.stats import chi2, poisson

from . import stats

FOLD_MAX = 1.5            # controls that disagree by more than this ...
CONTROL_P = 0.01          # ... and significantly, are flagged
P_MISS_MAX = 0.05
DISCARD_MAX = 0.5
DISCARD_IDENTITY_MIN = 85.0   # discarded reads this similar to the query may be true candidate reads
IDENTITY_GAP_MAX = 10.0
BREADTH_MIN = 0.6         # observed over expected breadth below this suggests a partial loss
MIN_READS_COVERAGE = 10
ASSIGNABLE_MIN = 0.2      # below this the candidate cannot be told from the rest of the family


def _mean(v):
    return sum(v) / len(v) if v else float("nan")


def control_agreement(counts, family):
    ctrl = family.controls
    n = [counts.by_gene.get(g, 0) for g in ctrl]
    L = [family.eff_length(g) for g in ctrl]
    if len(ctrl) < 2 or sum(n) == 0 or min(L) <= 0:
        return float("nan"), float("nan")
    rates = [x / l for x, l in zip(n, L)]
    fold = max(rates) / min(rates) if min(rates) > 0 else float("inf")
    tot, Lt = sum(n), sum(L)
    exp = [tot * l / Lt for l in L]
    x2 = sum((o - e) ** 2 / e for o, e in zip(n, exp))
    return fold, float(chi2.sf(x2, len(ctrl) - 1))


def leave_one_out(counts, family):
    ctrl = family.controls
    if len(ctrl) < 2:
        return float("nan"), float("nan")
    obs = counts.by_gene.get(family.candidate, 0)
    vals = []
    for g in ctrl:
        keep = [c for c in ctrl if c != g]
        if sum(family.eff_length(c) for c in keep) <= 0:
            continue
        e = stats.expected_reads([counts.by_gene.get(c, 0) for c in keep], [family.eff_length(c) for c in keep],
                                 family.eff_length(family.candidate))
        if e > 0:
            vals.append(obs / e)
    return (min(vals), max(vals)) if vals else (float("nan"), float("nan"))


def p_miss_one_copy(expected, th):
    """P(ratio <= absent threshold) when the true ratio is 0.5."""
    if not expected or expected <= 0:
        return float("nan")
    return float(poisson.cdf(math.floor(th.absent_max * expected), 0.5 * expected))


def coverage(spans, length):
    if not spans or not length:
        return float("nan"), float("nan"), float("nan")
    L = int(length)
    cov = [False] * (L + 2)
    for a, b in spans:
        for i in range(max(1, a), min(L, b) + 1):
            cov[i] = True
    breadth = sum(cov[1:L + 1]) / L
    mean_len = _mean([b - a + 1 for a, b in spans])
    expected = 1 - math.exp(-len(spans) * mean_len / L)
    gap = best = 0
    for i in range(1, L + 1):
        gap = 0 if cov[i] else gap + 1
        best = max(best, gap)
    return breadth, expected, best


def run(counts, family, res, th):
    fold, cp = control_agreement(counts, family)
    lo, hi = leave_one_out(counts, family)
    pm = p_miss_one_copy(res.expected, th)
    obs, disc = res.observed, counts.discarded
    dfrac = disc / (obs + disc) if obs + disc else float("nan")
    ctrl_id = _mean([x for g in family.controls for x in counts.identity.get(g, [])])
    gap = ctrl_id - _mean(counts.candidate_identity) if counts.candidate_identity else float("nan")
    br, bexp, lg = coverage(counts.candidate_spans, family.lengths[family.candidate])
    flags = []
    if fold == fold and fold > FOLD_MAX and cp < CONTROL_P:
        flags.append("controls_disagree")
    if lo == lo and res.status in (stats.ABSENT, stats.PRESENT) and \
            stats.classify(lo, res.expected, th) != stats.classify(hi, res.expected, th):
        flags.append("call_depends_on_one_control")
    if pm == pm and res.status == stats.ABSENT and pm > P_MISS_MAX:
        flags.append("low_power_for_one_copy")
    did = _mean(counts.discarded_identity)
    if dfrac == dfrac and disc >= 5 and dfrac > DISCARD_MAX and did >= DISCARD_IDENTITY_MIN:
        flags.append("rule_may_discard_true_reads")
    if gap == gap and obs >= 5 and gap > IDENTITY_GAP_MAX:
        flags.append("divergent_candidate")
    if br == br and obs >= MIN_READS_COVERAGE and bexp > 0 and br / bexp < BREADTH_MIN:
        flags.append("possible_partial_loss")
    ac = family.assignable.get(family.candidate)
    if ac is not None and ac < ASSIGNABLE_MIN:
        flags.append("low_assignable_fraction")
    return {"control_fold": fold, "control_p": cp, "ratio_loo_min": lo, "ratio_loo_max": hi,
            "p_miss_one_copy": pm, "discarded_fraction": dfrac, "discarded_identity": did, "identity_gap": gap,
            "breadth": br, "breadth_expected": bexp, "largest_gap_aa": lg,
            "assignable_candidate": ac if ac is not None else "", "flags": ",".join(flags)}
