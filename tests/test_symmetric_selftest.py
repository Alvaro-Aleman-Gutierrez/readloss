"""Symmetric rule (version 0.2.0) and selftest of detectability."""
import gzip
import os
import random
import shutil

import pytest

from readloss import assign, family, pipeline, selftest, stats

EX = os.path.join(os.path.dirname(__file__), "..", "ejemplos", "dcx")
KNOWN = {"Danionella_cerebrum_ERR3284980": stats.PRESENT, "Ictalurus_furcatus_SRR18966194": stats.PRESENT,
         "Paedocypris_micromegethes_SRR11971489": stats.ABSENT, "Paramisgurnus_dabryanus_SRR23104325": stats.ABSENT}


def two_gene_family(tmp_path, cand="cand", ctrl="ctrl"):
    (tmp_path / "g.tsv").write_text(f"gene\trole\tlength_aa\n{cand}\tcandidate\t100\n{ctrl}\tcontrol\t100\n")
    (tmp_path / "p.tsv").write_text(f"protein_id\tgene\npC\t{cand}\npK\t{ctrl}\n")
    return family.load(tmp_path / "g.tsv", tmp_path / "p.tsv")


def near_identical_hits(path, present=True, n=200, diagnostic=0.10, seed=1, prefix=""):
    """Reads of two copies that differ at few positions: most reads tie, a few differ by 5 points."""
    rng = random.Random(seed)
    with open(path, "w") as fh:
        for src in (("cand", "ctrl") if present else ("ctrl",)):
            for i in range(n):
                q = f"{prefix}{src}{selftest.SEP}{i}" if prefix is not None else f"{src}_{i}"
                diag = rng.random() < diagnostic
                sc_c = 150 if (src == "cand" or not diag) else 145
                sc_k = 150 if (src == "ctrl" or not diag) else 145
                fh.write(f"{q}\tpC\t99.5\t50\t{sc_c}\n{q}\tpK\t99.5\t50\t{sc_k}\n")
    return path


def run(path, fam_, mode):
    rule = assign.Rule(mode=mode)
    return assign.summarise(assign.assign(assign.read_hits(path, fam_, rule), fam_, rule), fam_, rule=rule)


def test_census_rule_calls_near_identical_copy_absent(tmp_path):
    f = two_gene_family(tmp_path)
    s = run(near_identical_hits(tmp_path / "h.tsv"), f, assign.CENSUS)
    assert s["status"] == stats.ABSENT          # the bias the symmetric rule removes


def test_symmetric_rule_near_identical_copy_present(tmp_path):
    f = two_gene_family(tmp_path)
    s = run(near_identical_hits(tmp_path / "h.tsv"), f, assign.SYMMETRIC)
    assert s["status"] in (stats.PRESENT, stats.UNINFORMATIVE)
    assert s["status"] != stats.ABSENT
    assert 0.5 < s["ratio"] < 2.0
    assert s["control_discarded"] > 100          # ties removed from the control side too


def test_symmetric_rule_near_identical_copy_absent(tmp_path):
    f = two_gene_family(tmp_path)
    s = run(near_identical_hits(tmp_path / "h.tsv", present=False, n=1000), f, assign.SYMMETRIC)
    assert s["observed"] == 0 and s["status"] == stats.ABSENT


def test_symmetric_rule_independent_of_gene_names(tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    a = two_gene_family(tmp_path / "a")
    b = two_gene_family(tmp_path / "b", cand="zzz", ctrl="aaa")
    h = near_identical_hits(tmp_path / "h.tsv")
    sa, sb = run(h, a, assign.SYMMETRIC), run(h, b, assign.SYMMETRIC)
    assert (sa["observed"], sa["expected"]) == (sb["observed"], sb["expected"])


def test_tie_between_candidate_and_other_never_counts(tmp_path):
    (tmp_path / "g.tsv").write_text("gene\trole\tlength_aa\ncand\tcandidate\t100\nctrl\tcontrol\t100\n"
                                    "para\tother\t0\n")
    (tmp_path / "p.tsv").write_text("protein_id\tgene\npC\tcand\npK\tctrl\npP\tpara\n")
    f = family.load(tmp_path / "g.tsv", tmp_path / "p.tsv")
    h = tmp_path / "h.tsv"
    h.write_text("r1\tpC\t95\t50\t100\nr1\tpP\t95\t50\t100\n"     # candidate = other: not the candidate
                 "r2\tpK\t95\t50\t100\nr2\tpP\t95\t50\t100\n"     # control = other: not the control
                 "r3\tpK\t95\t50\t100\nr3\tpP\t95\t50\t90\n")     # control wins: counts
    c = assign.assign(assign.read_hits(h, f), f, assign.Rule())
    assert c.by_gene["cand"] == 0 and c.by_gene["ctrl"] == 1 and c.by_gene["para"] == 2


def unzip(tmp, stem):
    p = tmp / f"{stem}.tsv"
    with gzip.open(os.path.join(EX, "hits", f"{stem}.hits.tsv.gz"), "rt") as a, open(p, "w") as b:
        shutil.copyfileobj(a, b)
    return p


@pytest.mark.parametrize("stem", sorted(KNOWN))
def test_symmetric_rule_keeps_doublecortin_calls(tmp_path, stem):
    f = family.load(os.path.join(EX, "genes.tsv"), os.path.join(EX, "proteins.tsv"))
    h = unzip(tmp_path, stem)
    sym, cen = run(h, f, assign.SYMMETRIC), run(h, f, assign.CENSUS)
    assert sym["status"] == cen["status"] == KNOWN[stem]
    assert abs(sym["ratio"] - cen["ratio"]) < 0.05


def test_fragments(tmp_path):
    f = two_gene_family(tmp_path)
    seqs = {"pC": "A" * 120, "pK": "C" * 60}
    fr = selftest.fragments(f, seqs, length=50, step=10)
    assert sum(1 for i, _ in fr if i.startswith("cand::")) == 8        # starts 1, 11, ..., 71
    assert sum(1 for i, _ in fr if i.startswith("ctrl::")) == 2
    assert all(len(s) == 50 for _, s in fr)
    mut = selftest.fragments(f, {"pC": "A" * 1000, "pK": "C" * 60}, length=1000, step=1000, divergence=0.2)
    changed = sum(1 for x in mut[0][1] if x != "A")
    assert 150 < changed < 250


def test_selftest_from_hits_and_correction(tmp_path):
    f = two_gene_family(tmp_path)
    sim = near_identical_hits(tmp_path / "sim.tsv", prefix="")
    rows = {r["gene"]: r for r in selftest.run_hits(f, sim, assign.Rule())}
    assert 0.05 < rows["cand"]["assignable_fraction"] < 0.2              # only diagnostic reads survive
    assert 0.05 < rows["ctrl"]["assignable_fraction"] < 0.2
    out = tmp_path / "self.tsv"
    out.write_text("gene\trole\tassignable_fraction\n" +
                   "".join(f"{g}\t{r['role']}\t{r['assignable_fraction']}\n" for g, r in rows.items()))
    family.load_assignable(out, f)
    s = run(near_identical_hits(tmp_path / "real.tsv", n=2000, seed=5), f, assign.SYMMETRIC)
    assert s["status"] == stats.UNINFORMATIVE and "low_assignable_fraction" in s["flags"]


def test_assignable_corrects_expected(tmp_path):
    f = two_gene_family(tmp_path)
    f.assignable = {"cand": 0.5, "ctrl": 1.0}          # candidate reads half as assignable as control reads
    h = tmp_path / "h.tsv"
    rows = [f"k{i}\tpK\t99\t50\t100\n" for i in range(100)] + [f"c{i}\tpC\t99\t50\t100\n" for i in range(50)]
    h.write_text("".join(rows))
    s = run(h, f, assign.SYMMETRIC)
    assert s["expected"] == pytest.approx(50.0) and s["ratio"] == pytest.approx(1.0)
    assert s["assignable_candidate"] == 0.5


def test_load_assignable_needs_all_genes(tmp_path):
    f = two_gene_family(tmp_path)
    t = tmp_path / "s.tsv"
    t.write_text("gene\trole\tassignable_fraction\ncand\tcandidate\t0.9\n")
    with pytest.raises(ValueError):
        family.load_assignable(t, f)


@pytest.mark.skipif(not (os.environ.get("READLOSS_DIAMOND") or shutil.which("diamond")),
                    reason="DIAMOND not installed")
def test_selftest_with_diamond(tmp_path):
    f = family.load(os.path.join(EX, "genes.tsv"), os.path.join(EX, "proteins.tsv"))
    rows = {r["gene"]: r for r in selftest.run_protein(f, os.path.join(EX, "family.faa"), assign.Rule(),
                                                       workdir=str(tmp_path), threads=2)}
    assert rows["DCX"]["assignable_fraction"] > 0.9
    assert rows["dclk1a"]["assignable_fraction"] > 0.8 and rows["dclk2a"]["assignable_fraction"] > 0.8
    # two identical proteins in different genes cannot be told apart
    seq = selftest.read_fasta(os.path.join(EX, "family.faa"))
    pid = next(p for p, g in f.protein_gene.items() if g == "DCX")
    fa = tmp_path / "twin.faa"
    fa.write_text(f">pC\n{seq[pid]}\n>pK\n{seq[pid]}\n")
    t = two_gene_family(tmp_path)
    twin = {r["gene"]: r for r in selftest.run_protein(t, str(fa), assign.Rule(), workdir=str(tmp_path / "w"),
                                                       threads=2)}
    assert twin["cand"]["assignable_fraction"] == 0.0
    assert pipeline.which("diamond", "READLOSS_DIAMOND")
