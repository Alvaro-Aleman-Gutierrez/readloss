"""Nucleotide mode, contamination flag and fallback for older DIAMOND (version 0.3.0)."""
import io
import os
import random
import shutil
import stat
from contextlib import redirect_stdout

import pytest

from readloss import assign, cli, family, pipeline, selftest, stats

HAS_BWA = bool(os.environ.get("READLOSS_BWA") or shutil.which("bwa-mem2"))


def two_genes(tmp_path, cand_id="hC", ctrl_id="hE"):
    (tmp_path / "g.tsv").write_text("gene\trole\tlength_aa\ncand\tcandidate\t500\nctrl\tcontrol\t500\n")
    (tmp_path / "p.tsv").write_text(f"protein_id\tgene\n{cand_id}\tcand\n{ctrl_id}\tctrl\n")
    return tmp_path / "g.tsv", tmp_path / "p.tsv"


# --- SAM conversion ----------------------------------------------------------------------------

def test_sam_to_hits():
    sam = ["@SQ\tSN:hC\tLN:1500\n",
           "r1\t0\thC\t31\t60\t150M\t*\t0\t0\t" + "A" * 150 + "\t*\tNM:i:2\tAS:i:140\n",
           "r1\t256\thE\t31\t0\t150M\t*\t0\t0\t*\t*\tNM:i:4\tAS:i:130\n",
           "r2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n",
           "r3\t16\thC\t1\t60\t10S100M2I38M\t*\t0\t0\t" + "G" * 150 + "\t*\tNM:i:2\tAS:i:120\n"]
    rows = list(pipeline.sam_to_hits(sam))
    assert len(rows) == 3
    q, s, pid, ln, bs, st, en, tr, seq = rows[0]
    assert (q, s, pid, ln, bs, st, en, tr) == ("r1", "hC", round(100 * 148 / 150, 2), 50, 140, 11, 60, "-")
    assert rows[1][8] == "*"                          # secondary: no sequence
    assert rows[2][3] == 140 // 3 and rows[2][5] == 1


# --- homeologs that encode the same protein ----------------------------------------------------

def homeologs(seed=7, n=1500, div=0.015):
    rng = random.Random(seed)
    c = "".join(rng.choice("ACGT") for _ in range(n))
    e = "".join(rng.choice("ACGT".replace(b, "")) if rng.random() < div else b for b in c)
    return c, e


@pytest.mark.skipif(not HAS_BWA, reason="bwa-mem2 not installed")
def test_homeologs_in_nucleotide_mode(tmp_path):
    g, p = two_genes(tmp_path)
    c, e = homeologs()
    cds = tmp_path / "cds.fna"
    cds.write_text(f">hC\n{c}\n>hE\n{e}\n")
    f = family.load(g, p)
    db = pipeline.makedb_nt(str(cds), str(tmp_path / "fam"))
    rule = assign.Rule(margin=5, min_identity=95)
    # selftest: most reads can be told apart at 1.5% divergence
    st = {r["gene"]: r for r in selftest.run_nt(f, str(cds), rule, step=5, workdir=str(tmp_path / "st"), threads=2)}
    assert st["cand"]["assignable_fraction"] > 0.6 and st["ctrl"]["assignable_fraction"] > 0.6
    seqs = {"hC": c, "hE": e}
    both = selftest.reads_nt(f, seqs, step=3, seed=11)                 # specimen with both copies
    only_e = [r for r in both if r[0].startswith("ctrl")]              # specimen that lost the C copy
    out = {}
    for name, rs in (("both", both), ("lost", only_e)):
        fa = tmp_path / f"{name}.fa"
        fa.write_text("".join(f">{i}\n{s}\n" for i, s in rs))
        tsv = pipeline.search_nt(str(fa), db, str(tmp_path / f"{name}.tsv"), threads=2)
        info = {}
        h = assign.read_hits(tsv, f, rule, info=info)
        out[name] = assign.summarise(assign.assign(h, f, rule, info=info), f, rule=rule)
    assert out["both"]["status"] == stats.PRESENT and 0.8 < out["both"]["ratio"] < 1.25
    assert out["lost"]["status"] == stats.ABSENT and out["lost"]["observed"] <= 2


@pytest.mark.skipif(not HAS_BWA, reason="bwa-mem2 not installed")
def test_cli_sam_input(tmp_path):
    g, p = two_genes(tmp_path)
    c, e = homeologs()
    sam = tmp_path / "a.sam"
    lines = []
    for i in range(60):
        lines.append(f"c{i}\t0\thC\t{1 + 20 * i}\t60\t150M\t*\t0\t0\t{c[20 * i:20 * i + 150]}\t*\tNM:i:0\tAS:i:150\n")
        lines.append(f"c{i}\t256\thE\t{1 + 20 * i}\t0\t150M\t*\t0\t0\t*\t*\tNM:i:2\tAS:i:140\n")
        lines.append(f"e{i}\t0\thE\t{1 + 20 * i}\t60\t150M\t*\t0\t0\t{e[20 * i:20 * i + 150]}\t*\tNM:i:0\tAS:i:150\n")
        lines.append(f"e{i}\t256\thC\t{1 + 20 * i}\t0\t150M\t*\t0\t0\t*\t*\tNM:i:2\tAS:i:140\n")
    sam.write_text("".join(lines))
    buf = io.StringIO()
    with redirect_stdout(buf):
        cli.main(["test", "--sam", str(sam), "--genes", str(g), "--proteins", str(p), "--workdir", str(tmp_path)])
    import csv
    r = list(csv.DictReader(io.StringIO(buf.getvalue()), delimiter="\t"))[0]
    assert r["status"] == "present" and float(r["ratio"]) == pytest.approx(1.0)


def test_protein_selftest_cannot_split_identical_proteins(tmp_path):
    # the reason for the nucleotide mode: the same protein, two genes
    g, p = two_genes(tmp_path)
    f = family.load(g, p)
    hits = tmp_path / "h.tsv"
    hits.write_text("".join(f"cand::{i}\thC\t100\t50\t100\ncand::{i}\thE\t100\t50\t100\n" for i in range(50)))
    rows = {r["gene"]: r for r in selftest.run_hits(f, hits, assign.Rule())}
    assert rows["cand"]["assignable_fraction"] == 0.0


# --- contamination -------------------------------------------------------------------------------

def fam3(tmp_path):
    (tmp_path / "g.tsv").write_text("gene\trole\tlength_aa\ncand\tcandidate\t100\na\tcontrol\t100\nb\tcontrol\t100\n")
    (tmp_path / "p.tsv").write_text("protein_id\tgene\nc1\tcand\na1\ta\nb1\tb\n")
    return family.load(tmp_path / "g.tsv", tmp_path / "p.tsv")


def run_rows(tmp_path, f, rows):
    p = tmp_path / "h.tsv"
    p.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))
    return assign.summarise(assign.assign(assign.read_hits(p, f), f), f)


def rows_of(prefix, prot, n, pid):
    return [[f"{prefix}{i}", prot, pid, 40, 80, 1 + (7 * i) % 60, 40 + (7 * i) % 60] for i in range(n)]


def test_contamination_flag(tmp_path):
    f = fam3(tmp_path)
    s = run_rows(tmp_path, f, rows_of("a", "a1", 80, 84) + rows_of("b", "b1", 80, 86) + rows_of("c", "c1", 3, 100))
    assert s["status"] == stats.ABSENT and "possible_contamination" in s["flags"]


def test_conserved_candidate_is_not_contamination(tmp_path):
    f = fam3(tmp_path)
    s = run_rows(tmp_path, f, rows_of("a", "a1", 80, 84) + rows_of("b", "b1", 80, 86) + rows_of("c", "c1", 80, 99))
    assert s["status"] == stats.PRESENT and "possible_contamination" not in s["flags"]


# --- older DIAMOND without the read-sequence columns -------------------------------------------------

def test_fallback_without_sequence_fields(tmp_path, monkeypatch):
    fake = tmp_path / "diamond"
    fake.write_text("#!/bin/sh\nfor x in \"$@\"; do [ \"$x\" = qseq_translated ] && exit 1; done\n"
                    "out=''; prev=''; for x in \"$@\"; do [ \"$prev\" = -o ] && out=$x; prev=$x; done\n"
                    "[ -n \"$out\" ] && echo \"$@\" > \"$out\"; exit 0\n")
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("READLOSS_DIAMOND", str(fake))
    pipeline._SEQ_FIELDS.clear()
    assert pipeline.seq_fields_supported(str(fake)) is False
    out = pipeline.search("r.fa", "db", str(tmp_path / "o.tsv"), threads=1)
    args = open(out).read()
    assert "send" in args and "full_qseq" not in args
    pipeline._SEQ_FIELDS.clear()
