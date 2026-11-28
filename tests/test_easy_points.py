"""Per-run copy numbers, PCR-free libraries, short reads and selftest in the sensitive mode (0.3.0)."""
import csv
import gzip
import io
import os
import shutil
from contextlib import redirect_stdout

import pytest

from readloss import assign, cli, family, pipeline

EX = os.path.join(os.path.dirname(__file__), "..", "ejemplos", "dcx")
D = os.path.join(os.path.dirname(__file__), "datos", "sensible")


def run(args):
    buf = io.StringIO()
    with redirect_stdout(buf):
        cli.main(args)
    return list(csv.DictReader(io.StringIO(buf.getvalue()), delimiter="\t"))


def test_parse_copies():
    assert cli.parse_copies("dclk1a=2,dclk2a=2") == {"dclk1a": 2.0, "dclk2a": 2.0}
    assert cli.parse_copies("a=2; b=1.5") == {"a": 2.0, "b": 1.5}
    with pytest.raises(ValueError):
        cli.parse_copies("a")
    with pytest.raises(ValueError):
        cli.parse_copies("a=0")


def test_copies_for_one_specimen(tmp_path):
    f = family.load(os.path.join(EX, "genes.tsv"), os.path.join(EX, "proteins.tsv"))
    g = cli._with_copies(f, "dclk1a=2,dclk2a=2")
    assert g.copies == {"dclk1a": 2.0, "dclk2a": 2.0} and f.copies == {}       # the shared family is untouched
    with pytest.raises(ValueError):
        cli._with_copies(f, "nosuchgene=2")


def test_copies_option_halves_expected(tmp_path):
    p = tmp_path / "h.tsv"
    with gzip.open(os.path.join(EX, "hits", "Danionella_cerebrum_ERR3284980.hits.tsv.gz"), "rt") as a, \
            open(p, "w") as b:
        shutil.copyfileobj(a, b)
    G = ["--genes", os.path.join(EX, "genes.tsv"), "--proteins", os.path.join(EX, "proteins.tsv")]
    one = run(["test", "--hits", str(p)] + G)[0]
    two = run(["test", "--hits", str(p), "--copies", "dclk1a=2,dclk2a=2"] + G)[0]
    assert float(two["expected"]) == pytest.approx(float(one["expected"]) / 2, abs=0.1)


def test_pcr_free():
    assert pipeline.pcr_free("TruSeq PCRfree") == "yes"
    assert pipeline.pcr_free("Illumina TruSeq DNA PCR-Free, 350 bp") == "yes"
    assert pipeline.pcr_free("amplification-free library") == "yes"
    assert pipeline.pcr_free("Nextera XT") == "no"
    assert pipeline.pcr_free("libraries were amplified with 8 cycles of PCR") == "no"
    assert pipeline.pcr_free("BDSW190028934") == "unknown"
    assert pipeline.pcr_free("") == "unknown"


def test_parse_experiment():
    xml = ("<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE><EXPERIMENT><DESIGN><DESIGN_DESCRIPTION>TruSeq PCRfree"
           "</DESIGN_DESCRIPTION><LIBRARY_DESCRIPTOR><LIBRARY_NAME>L1</LIBRARY_NAME></LIBRARY_DESCRIPTOR>"
           "</DESIGN></EXPERIMENT></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>")
    t = pipeline.parse_experiment(xml)
    assert "TruSeq PCRfree" in t and pipeline.pcr_free(t) == "yes"
    assert pipeline.parse_experiment("not xml") == ""


def test_metadata_reports_pcr_free():
    row = cli._apply_metadata({"status": "present", "flags": ""},
                              {"LibraryStrategy": "WGS", "LibrarySelection": "RANDOM", "LibrarySource": "GENOMIC",
                               "LibraryDescription": "TruSeq PCRfree"})
    assert row["pcr_free"] == "yes" and row["status"] == "present"


class A:
    min_aln = None


def test_short_reads_lower_min_aln():
    r = assign.Rule()
    assert cli._short_read_rule(r, A(), {"avgLength": "150", "LibraryLayout": "PAIRED"}).min_aln == 15  # 75 nt
    assert cli._short_read_rule(r, A(), {"avgLength": "90", "LibraryLayout": "SINGLE"}).min_aln == 18
    assert cli._short_read_rule(r, A(), {"avgLength": "200", "LibraryLayout": "PAIRED"}).min_aln == 30  # 100 nt
    assert cli._short_read_rule(r, A(), {"avgLength": "300", "LibraryLayout": "PAIRED"}).min_aln == 30
    given = A()
    given.min_aln = 30
    assert cli._short_read_rule(r, given, {"avgLength": "75", "LibraryLayout": "SINGLE"}).min_aln == 30
    assert cli._short_read_rule(r, A(), {}).min_aln == 30


def test_sensitive_accepts_assignable(tmp_path):
    bx = tmp_path / "b.tsv"
    with gzip.open(os.path.join(D, "cerebrum_ctrl.blastx.tsv.gz"), "rt") as a, open(bx, "w") as b:
        shutil.copyfileobj(a, b)
    st = tmp_path / "st.tsv"
    st.write_text("gene\trole\tassignable_fraction\nDCX\tcandidate\t0.5\ndclk1a\tcontrol\t1\n"
                  "dclk2a\tcontrol\t1\ndclk3\tcontrol\t1\n")
    G = ["--genes", os.path.join(EX, "genes_sensible.tsv"), "--proteins", os.path.join(EX, "proteins.tsv")]
    plain = run(["sensitive", "--blastx", str(bx)] + G)[0]
    corr = run(["sensitive", "--blastx", str(bx), "--assignable", str(st)] + G)[0]
    assert float(corr["expected"]) == pytest.approx(float(plain["expected"]) / 2, abs=0.1)


# --- found in the user pilot ---------------------------------------------------------------------

def test_write_keeps_header_columns(tmp_path, capsys):
    out = tmp_path / "r.tsv"
    cli._write([{"a": 1, "b": 2}], str(out))
    cli._write([{"a": 3, "b": 4, "c": 5}], str(out))
    rows = list(csv.DictReader(open(out), delimiter="\t"))
    assert [r["a"] for r in rows] == ["1", "3"] and "c" not in rows[0]
    assert "left out: c" in capsys.readouterr().err


def batch_args(tmp_path, runs):
    p = tmp_path / "runs.tsv"
    p.write_text(runs)
    return ["batch", "--runs", str(p), "--db", "x.dmnd", "--genes", os.path.join(EX, "genes.tsv"),
            "--proteins", os.path.join(EX, "proteins.tsv"), "--workdir", str(tmp_path / "w"),
            "--out", str(tmp_path / "res.tsv"), "--spots-max", "0"]


def test_batch_skips_non_wgs_and_uses_spots_per_run(tmp_path, monkeypatch):
    meta = {"RNA1": {"LibraryStrategy": "RNA-Seq", "LibrarySelection": "cDNA", "LibrarySource": "TRANSCRIPTOMIC"},
            "WGS1": {"LibraryStrategy": "WGS", "LibrarySelection": "RANDOM", "LibrarySource": "GENOMIC"}}
    monkeypatch.setattr(cli, "_metadata", lambda run, a: dict(meta[run]))
    asked = {}
    hits = tmp_path / "h.tsv"
    with gzip.open(os.path.join(EX, "hits", "Danionella_cerebrum_ERR3284980.hits.tsv.gz"), "rt") as a, \
            open(hits, "w") as b:
        shutil.copyfileobj(a, b)

    def fake_pass(run, spots, a, f, rule, th, wd, info=None, sex=None):
        asked[run] = spots
        row = cli._from_hits(str(hits), f, rule, th, a, sex=sex)
        row.update({"spots": spots, "hits_table": str(hits)})
        return row
    monkeypatch.setattr(cli, "_one_pass", fake_pass)
    cli.main(batch_args(tmp_path, "run\tlabel\tspots\nRNA1\trna\t12000000\nWGS1\twgs\t7000000\n"))
    rows = {r["run"]: r for r in csv.DictReader(open(tmp_path / "res.tsv"), delimiter="\t")}
    assert "RNA1" not in asked and asked["WGS1"] == 7000000          # no download for RNA-seq
    assert rows["RNA1"]["status"] == "uninformative" and "library_not_wgs" in rows["RNA1"]["flags"]
    assert rows["WGS1"]["status"] == "present" and rows["WGS1"]["library_strategy"] == "WGS"
    assert rows["RNA1"]["library_issues"] and rows["WGS1"]["library_issues"] == ""
