"""Diagnostics on synthetic hit tables."""
import math

import pytest

from readloss import assign, diagnostics, family, stats


@pytest.fixture
def fam(tmp_path):
    (tmp_path / "g.tsv").write_text("gene\trole\tlength_aa\ncand\tcandidate\t100\na\tcontrol\t100\n"
                                    "b\tcontrol\t100\npara\tother\t0\n")
    (tmp_path / "p.tsv").write_text("protein_id\tgene\nc1\tcand\na1\ta\nb1\tb\np1\tpara\n")
    return family.load(tmp_path / "g.tsv", tmp_path / "p.tsv")


def table(tmp_path, rows):
    p = tmp_path / "h.tsv"
    p.write_text("".join("\t".join(str(x) for x in r) + "\n" for r in rows))
    return p


def summary(tmp_path, fam, rows):
    return assign.summarise(assign.assign(assign.read_hits(table(tmp_path, rows), fam), fam), fam)


def reads(prefix, prot, n, pid=95, start=1, step=1, five=False):
    out = []
    for i in range(n):
        s = start + (i * step) % 60
        r = [f"{prefix}{i}", prot, pid, 40, 80]
        out.append(r if five else r + [s, s + 39])
    return out


def test_agreeing_controls_no_flags(tmp_path, fam):
    s = summary(tmp_path, fam, reads("a", "a1", 50) + reads("b", "b1", 50) + reads("c", "c1", 50, step=7))
    assert s["status"] == stats.PRESENT and s["flags"] == ""
    assert s["control_fold"] == pytest.approx(1.0)


def test_controls_disagree(tmp_path, fam):
    s = summary(tmp_path, fam, reads("a", "a1", 100) + reads("b", "b1", 30) + reads("c", "c1", 65, step=7))
    assert "controls_disagree" in s["flags"]


def test_partial_loss(tmp_path, fam):
    # 40 candidate reads all on the first 45 aa of a 100 aa protein
    s = summary(tmp_path, fam, reads("a", "a1", 40) + reads("b", "b1", 40) + reads("c", "c1", 40, step=0))
    assert s["breadth"] < 0.5 and "possible_partial_loss" in s["flags"]


def test_rule_may_discard_true_reads(tmp_path, fam):
    rows = reads("a", "a1", 50) + reads("b", "b1", 50)
    for i in range(20):        # high-identity candidate reads with a paralog hit within the margin
        rows += [[f"c{i}", "c1", 92, 40, 80, 1, 40], [f"c{i}", "p1", 90, 40, 78, 1, 40]]
    s = summary(tmp_path, fam, rows)
    assert s["candidate_discarded"] == 20 and "rule_may_discard_true_reads" in s["flags"]


def test_paralog_fragments_not_flagged(tmp_path, fam):
    rows = reads("a", "a1", 50) + reads("b", "b1", 50)
    for i in range(20):        # low-identity fragments: discarded, as intended, no flag
        rows += [[f"c{i}", "c1", 65, 40, 60, 1, 40], [f"c{i}", "p1", 64, 40, 58, 1, 40]]
    s = summary(tmp_path, fam, rows)
    assert s["status"] == stats.ABSENT and s["flags"] == ""


def test_five_columns_have_no_coverage(tmp_path, fam):
    s = summary(tmp_path, fam, reads("a", "a1", 50, five=True) + reads("b", "b1", 50, five=True)
                + reads("c", "c1", 50, five=True))
    assert math.isnan(s["breadth"]) and s["status"] == stats.PRESENT


def test_p_miss_one_copy():
    th = stats.Thresholds()
    assert diagnostics.p_miss_one_copy(10, th) == pytest.approx(0.0404, abs=1e-3)   # P(X <= 1 | 5)
    assert diagnostics.p_miss_one_copy(80, th) < 1e-8          # P(X <= 8 | 40), about 9e-10
