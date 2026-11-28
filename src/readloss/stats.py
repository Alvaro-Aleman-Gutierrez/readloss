"""Observed/expected read ratio, exact Poisson limits and classification.

The expected number of candidate reads is the control read rate per amino acid
times the candidate length:

    expected = (sum of control reads / sum of control lengths in aa) * candidate length

Confidence limits for the observed count use the exact (Garwood, 1936) chi-square
form and are divided by the expected count, which is treated as fixed.
"""
from dataclasses import dataclass, asdict
import math

from scipy.stats import chi2, poisson

ABSENT = "absent"
PRESENT = "present"
INTERMEDIATE = "intermediate"
UNINFORMATIVE = "uninformative"


@dataclass
class Thresholds:
    absent_max: float = 0.10       # ratio <= this is scored as absence
    present_min: float = 0.30      # ratio >= this is scored as presence
    min_expected: float = 10.0     # fewer expected reads than this is uninformative


@dataclass
class Result:
    observed: int
    expected: float
    ratio: float
    ci_low: float
    ci_high: float
    p_le_obs_if_present: float
    status: str

    def as_dict(self):
        return asdict(self)


def expected_reads(control_counts, control_lengths, candidate_length):
    """Expected candidate reads from control counts and protein lengths (aa)."""
    total_len = sum(control_lengths)
    if total_len <= 0:
        raise ValueError("control lengths must sum to a positive number")
    return sum(control_counts) / total_len * candidate_length


def garwood(k, alpha=0.05):
    """Exact two-sided Poisson confidence limits for a count k."""
    lo = chi2.ppf(alpha / 2, 2 * k) / 2 if k > 0 else 0.0
    hi = chi2.ppf(1 - alpha / 2, 2 * (k + 1)) / 2
    return lo, hi


def classify(ratio, expected, th=Thresholds()):
    if not expected or expected < th.min_expected or math.isnan(ratio):
        return UNINFORMATIVE
    if ratio <= th.absent_max:
        return ABSENT
    if ratio >= th.present_min:
        return PRESENT
    return INTERMEDIATE


def evaluate(observed, expected, th=Thresholds(), alpha=0.05):
    """Ratio, limits, probability of the observed count under presence and status."""
    if expected > 0:
        ratio = observed / expected
        lo, hi = garwood(observed, alpha)
        lo, hi = lo / expected, hi / expected
        p = float(poisson.cdf(observed, expected))
    else:
        ratio = lo = hi = p = float("nan")
    return Result(observed, expected, ratio, lo, hi, p, classify(ratio, expected, th))


# --- comparison of two Poisson rates, conditional on the total (version 0.2.0) ---
#
# The limits above treat the expected count as fixed, but it is estimated from the
# control reads, which are random too; with 30 control reads its coefficient of
# variation is about 18%. Conditioning on the total n = x + y of candidate and control
# reads, the candidate count is binomial with p = R*Lc / (R*Lc + Lk), where R is the
# ratio and Lc and Lk the effective lengths. Exact (Clopper-Pearson) limits for p give
# exact limits for R that include the uncertainty of both counts.

from scipy.stats import beta as _beta, binom as _binom

CN_STATES = (0.0, 0.5, 1.0, 2.0)
CN_FLOOR = 0.01          # ratio standing for "absent": reads from errors, contamination or paralogs


def _p_of_ratio(r, lc, lk):
    return r * lc / (r * lc + lk)


def _ratio_of_p(p, lc, lk):
    if p >= 1:
        return float("inf")
    return p / (1 - p) * lk / lc


def conditional(x, y, lc, lk, alpha=0.05):
    """Exact limits of the ratio and P(X <= x) under ratio 1, conditional on x + y.

    x, y: candidate and summed control reads; lc, lk: effective lengths.
    """
    n = x + y
    if n == 0 or lc <= 0 or lk <= 0:
        return float("nan"), float("nan"), float("nan")
    plo = _beta.ppf(alpha / 2, x, n - x + 1) if x > 0 else 0.0
    phi = _beta.ppf(1 - alpha / 2, x + 1, n - x) if x < n else 1.0
    p1 = _p_of_ratio(1.0, lc, lk)
    return _ratio_of_p(plo, lc, lk), _ratio_of_p(phi, lc, lk), float(_binom.cdf(x, n, p1))


def copy_number_states(x, y, lc, lk, states=CN_STATES, floor=CN_FLOOR):
    """Log-likelihood of each ratio state (0 stands for ``floor``), conditional on x + y.

    Returns (best state, log10 likelihood ratio of the best over the second best,
    {state: log10 likelihood}).
    """
    n = x + y
    if n == 0 or lc <= 0 or lk <= 0:
        return float("nan"), float("nan"), {}
    ll = {}
    for s in states:
        p = _p_of_ratio(max(s, floor), lc, lk)
        ll[s] = float(_binom.logpmf(x, n, p)) / math.log(10)
    ranked = sorted(ll, key=ll.get, reverse=True)
    return ranked[0], ll[ranked[0]] - ll[ranked[1]], ll
