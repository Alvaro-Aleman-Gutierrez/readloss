"""Reference coordinates and the figures of readloss plot."""
import csv
import gzip
import io
import os
import shutil
from contextlib import redirect_stdout

import pytest

from readloss import assign, cli, family

EX = os.path.join(os.path.dirname(__file__), "..", "ejemplos", "dcx")
G = ["--genes", os.path.join(EX, "genes.tsv"), "--proteins", os.path.join(EX, "proteins.tsv")]
pytest.importorskip("matplotlib")


def fam_with_coords():
    f = family.load(os.path.join(EX, "genes.tsv"), os.path.join(EX, "proteins.tsv"))
    return family.load_coords(os.path.join(EX, "family.faa"), f)


def test_reference_is_closest_length_and_maps_to_itself():
    f = fam_with_coords()
    ref = f.reference["DCX"]
    assert ref == "Caraur_LOC113084161_113084161"            # 357 aa, the length in genes.tsv
    assert [f.to_reference(ref, p) for p in (1, 100, 357)] == [1, 100, 357]


def test_longer_ortholog_maps_inside_the_reference():
    f = fam_with_coords()
    assert f.to_reference("Homsap_DCX_1641", 441) == 357     # human isoform with an N-terminal extension
    assert all(1 <= f.to_reference("Homsap_DCX_1641", p) <= 357 for p in range(1, 442, 20))


def test_breadth_is_bounded_with_coordinates(tmp_path):
    p = tmp_path / "h.tsv"
    with gzip.open(os.path.join(EX, "hits", "Ictalurus_furcatus_SRR18966194.hits.tsv.gz"), "rt") as a, \
            open(p, "w") as b:
        shutil.copyfileobj(a, b)
    f = fam_with_coords()
    c = assign.assign(assign.read_hits(p, f), f)
    assert all(1 <= a_ <= b_ <= 357 for a_, b_ in c.spans["DCX"])
    assert c.margins and all(g in ("DCX", "dclk1a", "dclk2a") for g, *_ in c.margins)


def run_plot(args):
    buf = io.StringIO()
    with redirect_stdout(buf):
        cli.main(["plot"] + args)
    return buf.getvalue().split()


@pytest.mark.parametrize("kind", ["coverage", "margins", "identity"])
def test_per_run_plots(tmp_path, kind):
    p = tmp_path / "h.tsv"
    with gzip.open(os.path.join(EX, "hits", "Paedocypris_micromegethes_SRR11971489.hits.tsv.gz"), "rt") as a, \
            open(p, "w") as b:
        shutil.copyfileobj(a, b)
    out = tmp_path / f"{kind}.png"
    png, tsv = run_plot([kind, "--hits", str(p), "--family-fasta", os.path.join(EX, "family.faa"),
                         "--out", str(out)] + G)
    assert os.path.getsize(png) > 10000
    rows = list(csv.DictReader(open(tsv), delimiter="\t"))
    assert rows
    if kind == "coverage":
        dcx = [r for r in rows if r["gene"] == "DCX"]
        assert len(dcx) == 357 and max(int(r["counted"]) for r in dcx) <= 1   # the gene is absent here


def test_selftest_and_power_plots(tmp_path):
    st = tmp_path / "st.tsv"
    st.write_text("gene\trole\tfragments\twith_hits\tassigned_self\tassigned_candidate\tassigned_other_controls\t"
                  "assigned_other_genes\tdiscarded\tno_hit\tassignable_fraction\n"
                  "DCX\tcandidate\t100\t100\t98\t0\t0\t0\t2\t0\t0.98\n"
                  "dclk1a\tcontrol\t100\t100\t90\t0\t0\t8\t2\t0\t0.90\n")
    png, tsv = run_plot(["selftest", "--selftest", str(st), "--out", str(tmp_path / "s.png")])
    assert os.path.getsize(png) > 5000 and len(list(csv.DictReader(open(tsv), delimiter="\t"))) == 12
    res = tmp_path / "r.tsv"
    res.write_text("label\texpected\tratio\tci_low_cond\tci_high_cond\tstatus\tflags\n"
                   "a\t16\t0.9\t0.5\t1.6\tpresent\t\nb\t4\t0.0\t0\t0.9\tuninformative\t\n"
                   "c\t18\t0.0\t0\t0.2\tabsent\t\nd\t0\tnan\tnan\tnan\tuninformative\tlibrary_not_wgs\n")
    png, tsv = run_plot(["power", "--results", str(res), "--out", str(tmp_path / "p.png")])
    assert len(list(csv.DictReader(open(tsv), delimiter="\t"))) == 3      # the run without expected reads is left out


def test_per_run_plot_needs_a_hit_table(tmp_path):
    res = tmp_path / "r.tsv"
    res.write_text("label\trun\thits_table\tstatus\tflags\nrna\tSRR1\t\tuninformative\tlibrary_not_wgs\n")
    with pytest.raises(SystemExit):
        run_plot(["coverage", "--results", str(res), "--run", "SRR1", "--out", str(tmp_path / "c.png")] + G)
