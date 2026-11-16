"""The statistics reproduce Supplementary Table S2 of the doublecortin study."""
import csv
import os

import pytest

from readloss import stats

D = os.path.join(os.path.dirname(__file__), "datos")
L = {"DCX": 357, "dclk1a": 770, "dclk2a": 826, "dclk3": 721}


def rows(name):
    with open(os.path.join(D, name)) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.mark.parametrize("r", rows("Table_S2.tsv"), ids=lambda r: r["species"])
def test_table_s2(r):
    ctrl = [int(r["reads_dclk1a"]), int(r["reads_dclk2a"]), int(r["reads_dclk3"])]
    exp = stats.expected_reads(ctrl, [L["dclk1a"], L["dclk2a"], L["dclk3"]], L["DCX"])
    res = stats.evaluate(int(r["reads_DCX"]), exp)
    assert round(exp, 1) == pytest.approx(float(r["expected_DCX"]))
    assert round(res.ratio, 3) == pytest.approx(float(r["obs_exp"]))
    assert round(res.ci_low, 3) == pytest.approx(float(r["CI95_low"]))
    assert round(res.ci_high, 3) == pytest.approx(float(r["CI95_high"]))
    assert float(f"{res.p_le_obs_if_present:.2g}") == pytest.approx(float(r["P_le_obs_if_present"]))
    assert res.status == r["DCX_status"]


def test_garwood_zero():
    lo, hi = stats.garwood(0)
    assert lo == 0 and hi == pytest.approx(3.689, abs=1e-3)


def test_uninformative():
    assert stats.evaluate(0, 5.0).status == stats.UNINFORMATIVE
    assert stats.evaluate(0, 0.0).status == stats.UNINFORMATIVE
