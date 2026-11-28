"""Diagnostics added to every result, so that users can judge a call before trusting it.

control_fold         highest over lowest reads per aa among the controls (1 = perfect agreement)
control_p            chi-square test that the controls collect reads in proportion to their length
ratio_loo_min/max    ratio recomputed leaving out each control in turn
p_miss_one_copy      probability of a false absence call if the candidate were present in one copy
                     while the controls have two (ratio 0.5), given the expected reads
discarded_fraction   candidate reads discarded by the margin or identity rule, over all candidate reads
discarded_identity   mean identity of the discarded reads (low = paralog fragments, high = possibly true reads)
identity_only_discarded
                     candidate reads that beat every known member but were discarded for low identity:
                     fragments of a gene missing from the family (flagged only if kept reads look the same)
identity_gap         mean identity of control reads minus mean identity of candidate reads
breadth              fraction of the candidate protein covered by kept reads (needs sstart and send)
breadth_expected     breadth expected for that number of reads if they fell at random along the protein
largest_gap_aa       longest stretch of the candidate protein without reads
assignable_candidate fraction of candidate reads the rule can assign, from selftest (blank without it)
expected_if_absent   candidate reads expected if the candidate were missing: control reads that the rule
                     would give to the candidate, at the rate measured by selftest (blank without it)
p_leakage            probability of at least the observed candidate reads from that leakage alone
stop_fraction_candidate / stop_fraction_controls
                     counted reads whose translated alignment holds a stop codon (needs 0.2.0 hit tables)
gc_candidate, gc_controls
                     GC of the coding sequences (--cds) or, without them, of the counted reads
gc_factor            coverage expected at the GC of the candidate over that at the GC of the controls,
                     from the panel (1 = no bias)
expected_panel       expected candidate reads from the median rate of the panel genes
control_vs_panel     control reads per aa over the panel median (about 1 for single-copy controls)
mean_read_length     mean length (nt) of the counted reads
sex, sex_source      sex used for sex-linked genes and where it came from (given, metadata, markers)
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
ASSIGNABLE_MIN = 0.5      # below this the candidate cannot be told from the rest of the family. Set from real
                          # data: X-Y gametologs of placental mammals at 79% protein identity (UTY, selftest 0.71)
                          # gave 13 of 14 correct calls; at 88-96% (DDX3Y 0.37, ZFY 0.34) females looked present
                          # and males looked absent in protein mode, while nucleotide mode called both correctly.
NT_RULE_IDENTITY = 90.0   # a minimum identity this high means nucleotide mode
NT_OWN_SPECIES = 98.5     # nucleotide mode: reads of a gene whose coding sequence comes from the same species
                          # match it at 99-100% (sequencing errors only). Male horse reads of Y genes, queried
                          # with cattle, pig, sheep, dog and human sequences, matched at 95.6-98.0% and half were
                          # lost to the 95% threshold, which gave false intermediate ratios.
LEAK_P = 0.05             # candidate reads this likely under the leakage of the controls alone: no presence call
STOP_MIN_FRACTION = 0.10  # candidate reads with stop codons above this ...
STOP_EXCESS = 3.0         # ... and this many times the control fraction suggest a pseudogene
# Reads discarded only for low identity, although they beat every known member by the margin, come from
# a source missing from the family; bunched in one region (breadth well below random placement) they
# point to a shared domain of an unlisted paralog.
GC_DIFF_MAX = 0.10        # GC difference between candidate and controls (Illumina PCR bias rises fast
GC_HIGH = 0.60            # above about 60% GC)
GC_FACTOR_RANGE = (0.7, 1.4)
PANEL_RANGE = (0.6, 1.6)  # control rate over panel rate outside this: controls not single copy
SHORT_READ_NT = 100       # reads shorter than this leave few alignments of 30 aa
# Contamination: reads of a gene present in the specimen diverge from the queries as much as the control
# reads, because both come from the same genome. A few candidate reads far closer to the queries than the
# control reads (and nearly identical to them) more likely come from another, closer species: laboratory
# contamination, index hopping between samples of one sequencing run, or DNA of the reference species.
# A strongly conserved candidate also gives closer reads, but then the ratio is not low, so the flag is
# raised only when the ratio is below the presence threshold.
CONTAM_GAP = 10.0         # candidate identity minus control identity, percentage points
CONTAM_IDENTITY = 97.0    # mean identity of the candidate reads to the queries


def _mean(v):
    return sum(v) / len(v) if v else float("nan")


def control_agreement(counts, family, sex="unknown"):
    ctrl = family.controls
    n = [counts.by_gene.get(g, 0) for g in ctrl]
    L = [family.eff_length(g, sex) for g in ctrl]
    if len(ctrl) < 2 or sum(n) == 0 or min(L) <= 0:
        return float("nan"), float("nan")
    rates = [x / l for x, l in zip(n, L)]
    fold = max(rates) / min(rates) if min(rates) > 0 else float("inf")
    tot, Lt = sum(n), sum(L)
    exp = [tot * l / Lt for l in L]
    x2 = sum((o - e) ** 2 / e for o, e in zip(n, exp))
    return fold, float(chi2.sf(x2, len(ctrl) - 1))


def leave_one_out(counts, family, sex="unknown"):
    ctrl = family.controls
    if len(ctrl) < 2:
        return float("nan"), float("nan")
    obs = counts.by_gene.get(family.candidate, 0)
    vals = []
    for g in ctrl:
        keep = [c for c in ctrl if c != g]
        if sum(family.eff_length(c, sex) for c in keep) <= 0:
            continue
        e = stats.expected_reads([counts.by_gene.get(c, 0) for c in keep], [family.eff_length(c, sex) for c in keep],
                                 family.eff_length(family.candidate, sex))
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


def _frac(a, b):
    return a / b if b else float("nan")


def run(counts, family, res, th, ctx=None):
    ctx = ctx or {}
    sex = ctx.get("sex", "unknown")
    fold, cp = control_agreement(counts, family, sex)
    lo, hi = leave_one_out(counts, family, sex)
    pm = p_miss_one_copy(res.expected, th)
    obs, disc = res.observed, counts.discarded
    dfrac = disc / (obs + disc) if obs + disc else float("nan")
    # control reads are not filtered by identity, so reads of distant members of a large family (other
    # DEAD-box helicases, zinc fingers, kinases) pull their mean down; compare with those that would pass
    # the candidate threshold
    mi = ctx.get("min_identity", 0.0) or 0.0
    ctrl_id = _mean([x for g in family.controls for x in counts.identity.get(g, []) if x >= mi])
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
    cand = family.candidate
    # pseudogene: stop codons inside the alignments of candidate reads, beyond the sequencing-error background
    sc = _frac(counts.stops.get(cand, 0), counts.with_seq.get(cand, 0))
    nk = sum(counts.with_seq.get(g, 0) for g in family.controls)
    sk = _frac(sum(counts.stops.get(g, 0) for g in family.controls), nk)
    if sc == sc and counts.with_seq.get(cand, 0) >= 5 and sc >= STOP_MIN_FRACTION and \
            sc > STOP_EXCESS * (sk if sk == sk else 0) + 0.02:
        flags.append("possible_pseudogene")
    # a family member missing from the database: kept candidate reads less similar than the controls and
    # bunched in one region, or many discarded reads of low identity bunched in one region
    bunched = br == br and bexp > 0 and br / bexp < BREADTH_MIN
    io = counts.identity_only_spans
    dbr, dbexp, _ = coverage(io, family.lengths[cand])
    disc_bunched = dbr == dbr and dbexp > 0 and dbr / dbexp < BREADTH_MIN
    # only a threat when kept candidate reads could be such fragments too: they are bunched in one region and
    # either less similar than the controls or accompanied by bunched identity-only discards. Discarded
    # fragments alone (as in a species that lacks the gene) leave the call safe and are only reported.
    if obs >= 5 and bunched and ((gap == gap and gap > IDENTITY_GAP_MAX) or (len(io) >= 5 and disc_bunched)):
        flags.append("possible_missing_paralog")
    gc_c, gc_k, gf = ctx.get("gc_candidate", float("nan")), ctx.get("gc_controls", float("nan")), \
        ctx.get("gc_factor", float("nan"))
    if gc_c == gc_c and gc_k == gc_k and (abs(gc_c - gc_k) > GC_DIFF_MAX or
                                          (max(gc_c, gc_k) > GC_HIGH and abs(gc_c - gc_k) > GC_DIFF_MAX / 2)):
        flags.append("gc_mismatch")
    if gf == gf and not GC_FACTOR_RANGE[0] <= gf <= GC_FACTOR_RANGE[1]:
        flags.append("gc_bias")
    prate = ctx.get("panel_rate", float("nan"))
    ctrl_n = sum(counts.by_gene.get(g, 0) for g in family.controls)
    ctrl_l = sum(family.eff_length(g, sex) for g in family.controls)
    cvp = (ctrl_n / ctrl_l) / prate if prate == prate and prate > 0 and ctrl_l > 0 else float("nan")
    exp_panel = prate * ctx.get("lc", 0) if prate == prate else float("nan")
    if cvp == cvp and not PANEL_RANGE[0] <= cvp <= PANEL_RANGE[1]:
        flags.append("controls_copy_number")
    if exp_panel == exp_panel and exp_panel >= th.min_expected and res.status in (stats.ABSENT, stats.PRESENT):
        if stats.classify(obs / exp_panel, exp_panel, th) != res.status:
            flags.append("panel_disagrees")
    cid = _mean(counts.candidate_identity)
    if obs >= 2 and cid == cid and ctrl_id == ctrl_id and cid - ctrl_id > CONTAM_GAP and \
            cid >= CONTAM_IDENTITY and res.ratio == res.ratio and res.ratio < th.present_min:
        flags.append("possible_contamination")
    if mi >= NT_RULE_IDENTITY and obs >= 3 and cid == cid and cid < NT_OWN_SPECIES:
        flags.append("candidate_reference_divergent")
    rl = _mean(counts.read_len)
    if rl == rl and rl < SHORT_READ_NT:
        flags.append("short_reads")
    linked = family.sex_linked(cand) or any(family.sex_linked(g) for g in family.controls)
    if linked and sex == "unknown":
        flags.append("sex_unknown_for_sex_linked_gene")
    if ctx.get("lc", 1) <= 0:
        flags.append("not_expected_in_this_sex")
    return {"control_fold": fold, "control_p": cp, "ratio_loo_min": lo, "ratio_loo_max": hi,
            "p_miss_one_copy": pm, "discarded_fraction": dfrac, "discarded_identity": did, "identity_gap": gap,
            "breadth": br, "breadth_expected": bexp, "largest_gap_aa": lg,
            "assignable_candidate": ac if ac is not None else "",
            "identity_only_discarded": len(io), "stop_fraction_candidate": sc, "stop_fraction_controls": sk,
            "gc_candidate": gc_c, "gc_controls": gc_k, "gc_factor": gf,
            "expected_panel": exp_panel, "control_vs_panel": cvp, "panel_genes": ctx.get("panel_genes", 0),
            "mean_read_length": rl, "sex": sex, "sex_source": ctx.get("sex_source", ""),
            "marker_ratio": ctx.get("marker_ratio", float("nan")),
            "flags": ",".join(flags)}
