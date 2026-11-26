"""Command line on the example hit tables (no download)."""
import csv
import gzip
import io
import os
import shutil
from contextlib import redirect_stdout

from readloss import cli

EX = os.path.join(os.path.dirname(__file__), "..", "ejemplos", "dcx")
G = ["--genes", os.path.join(EX, "genes.tsv"), "--proteins", os.path.join(EX, "proteins.tsv"), "--rule", "census"]


def unzip(tmp, stem):
    p = tmp / f"{stem}.tsv"
    with gzip.open(os.path.join(EX, "hits", f"{stem}.hits.tsv.gz"), "rt") as a, open(p, "w") as b:
        shutil.copyfileobj(a, b)
    return str(p)


def run(args):
    buf = io.StringIO()
    with redirect_stdout(buf):
        cli.main(args)
    return list(csv.DictReader(io.StringIO(buf.getvalue()), delimiter="\t"))


def test_test_hits(tmp_path):
    r = run(["test", "--hits", unzip(tmp_path, "Paedocypris_micromegethes_SRR11971489")] + G)[0]
    assert r["status"] == "absent" and r["observed"] == "1" and r["expected"] == "74.0"


def test_depth(tmp_path):
    r = run(["depth", "--hits", unzip(tmp_path, "Paramisgurnus_dabryanus_SRR23104325"),
             "--spots", "35000000"] + G)[0]
    assert abs(float(r["expected_per_million_spots"]) - 0.537) < 0.01
    assert r["recommended_spots"] == "28000000"          # 1.5 x 10 / 0.537 per M, rounded up to whole millions


def test_calibrate(tmp_path):
    s = tmp_path / "s.tsv"
    s.write_text("label\thits\tknown_status\n"
                 f"a\t{unzip(tmp_path, 'Danionella_cerebrum_ERR3284980')}\tpresent\n"
                 f"b\t{unzip(tmp_path, 'Paedocypris_micromegethes_SRR11971489')}\tabsent\n")
    rows = run(["calibrate", "--samples", str(s), "--margins", "5", "--identities", "75"] + G)
    sep = [r for r in rows if r["label"] == "SEPARATION"][0]
    assert sep["status"] == "separated" and abs(float(sep["ratio"]) - 0.774) < 1e-3
