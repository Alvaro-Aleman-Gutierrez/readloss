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
