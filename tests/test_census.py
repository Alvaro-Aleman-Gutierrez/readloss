"""Classification reproduces the read census (Supplementary Table S3) of the doublecortin study."""
import csv
import os

import pytest

from readloss import stats

D = os.path.join(os.path.dirname(__file__), "datos")
L = {"DCX": 357, "dclk1a": 770, "dclk2a": 826}
MAP = {"present": stats.PRESENT, "absent": stats.ABSENT, "ambiguous": stats.INTERMEDIATE}


def rows(name):
    with open(os.path.join(D, name)) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.mark.parametrize("r", rows("census_read_counts_per_run.tsv"), ids=lambda r: r["run"])
def test_counts_to_ratio(r):
    exp = stats.expected_reads([int(r["reads_dclk1a"]), int(r["reads_dclk2a"])],
                               [L["dclk1a"], L["dclk2a"]], L["DCX"])
    assert round(exp, 1) == pytest.approx(float(r["expected_DCX"]))
    assert round(int(r["reads_DCX"]) / exp, 3) == pytest.approx(float(r["obs_exp"]), abs=1e-3)


@pytest.mark.parametrize("r", rows("Table_S3.tsv"), ids=lambda r: r["run"])
def test_status(r):
    st = stats.classify(float(r["obs_exp"]), float(r["expected_DCX"]))
    assert st == MAP[r["DCX_status"]]
    if r["obs_exp_15M"]:            # repeated at 40 M spots, so the first pass must have been intermediate
        assert stats.classify(float(r["obs_exp_15M"]), float(r["expected_DCX_15M"])) == stats.INTERMEDIATE
