"""Diagnostics and corrections with a biological basis (version 0.2.0)."""
import math
import os
import shutil
import stat
import subprocess

import pytest

from readloss import assign, cli, extras, family, pipeline, selftest, stats

HAS_DIAMOND = bool(os.environ.get("READLOSS_DIAMOND") or shutil.which("diamond"))


def fam_files(tmp_path, rows, prots):
    (tmp_path / "g.tsv").write_text("gene\trole\tlength_aa\tcopies\tlocation\n" +
                                    "".join("\t".join(map(str, r)) + "\n" for r in rows))
    (tmp_path / "p.tsv").write_text("protein_id\tgene\n" + "".join(f"{p}\t{g}\n" for p, g in prots))
    return tmp_path / "g.tsv", tmp_path / "p.tsv"


def basic(tmp_path, cand_loc="", cand_copies="", ctrl_copies=""):
    g, p = fam_files(tmp_path, [("cand", "candidate", 100, cand_copies, cand_loc),
                                ("a", "control", 100, ctrl_copies, ""), ("b", "control", 100, ctrl_copies, ""),
                                ("para", "other", 0, "", "")],
                     [("c1", "cand"), ("a1", "a"), ("b1", "b"), ("p1", "para")])
    return family.load(g, p)


def write(tmp_path, rows, name="h.tsv"):
    p = tmp_path / name
    p.write_text("".join("\t".join(str(x) for x in r) + "\n" for r in rows))
    return p


def reads(prefix, prot, n, pid=95, bs=80, stop=False, seq=None, start=1, step=7, span=40):
    out = []
    for i in range(n):
        s = start + (i * step) % 60
        tr = "MKV*LLA" if stop else "MKVQLLA"
        out.append([f"{prefix}{i}", prot, pid, 40, bs, s, s + span - 1, tr, seq or ("ACGT" * 37 + "AC")])
    return out


def summary(tmp_path, f, rows, **kw):
    info = {}
    h = assign.read_hits(write(tmp_path, rows), f, info=info)
    return assign.summarise(assign.assign(h, f, info=info), f, **kw)


# --- SRA metadata ---------------------------------------------------------------------------

XML = """<?xml version="1.0" encoding="UTF-8"  ?>
<SraRunInfo><Row><Run>SRR1</Run><spots>1000</spots><avgLength>300</avgLength>
<LibraryStrategy>RNA-Seq</LibraryStrategy><LibrarySelection>cDNA</LibrarySelection>
<LibrarySource>TRANSCRIPTOMIC</LibrarySource><LibraryLayout>PAIRED</LibraryLayout>
<Sex>pooled male and female</Sex></Row></SraRunInfo>"""


def test_runinfo_xml_and_csv():
    i = pipeline.parse_runinfo(XML, "SRR1")
    assert i["LibraryStrategy"] == "RNA-Seq" and pipeline.mate_length(i) == 150
    assert pipeline.library_issues(i) == ["strategy RNA-SEQ", "selection CDNA", "source TRANSCRIPTOMIC"]
    c = pipeline.parse_runinfo("Run,LibraryStrategy,LibrarySelection,LibrarySource\nSRR2,WGS,RANDOM,GENOMIC\n",
                               "SRR2")
    assert pipeline.library_issues(c) == []
    assert pipeline.parse_runinfo(XML, "SRR9") == {}


def test_metadata_marks_rna_seq_and_pools():
    row = cli._apply_metadata({"status": "absent", "flags": ""}, pipeline.parse_runinfo(XML, "SRR1"))
    assert row["status"] == stats.UNINFORMATIVE
    assert "library_not_wgs" in row["flags"] and "pooled_sample" in row["flags"]


def test_capture_library_is_not_wgs():
    i = {"LibraryStrategy": "Targeted-Capture", "LibrarySelection": "Hybrid Selection", "LibrarySource": "GENOMIC"}
    assert len(pipeline.library_issues(i)) == 2


def test_spot_ranges():
    r = pipeline.spot_ranges(15_000_000, 60_000_000, 4)
    assert len(r) == 4 and sum(b - a + 1 for a, b in r) == 15_000_000
    assert r[1][0] == 15_000_001 and r[-1][1] <= 60_000_000
    assert pipeline.spot_ranges(15, 10, 4) == [(1, 15)]


def fake_dump(tmp_path):
    """A fastq-dump stand-in that writes one read per requested spot and reports them."""
    s = tmp_path / "fastq-dump"
    s.write_text("#!/bin/sh\nN=1; X=0\nwhile [ $# -gt 0 ]; do case $1 in -N) N=$2; shift;; -X) X=$2; shift;; esac; "
                 "shift; done\ni=$N; while [ $i -le $X ]; do echo \">r$i\"; echo ACGT; i=$((i+1)); done\n"
                 "echo \"Written $((X-N+1)) spots\" 1>&2\n")
    s.chmod(s.stat().st_mode | stat.S_IEXEC)
    return str(s)


def test_fetch_in_blocks(tmp_path, monkeypatch):
    monkeypatch.setenv("READLOSS_FASTQ_DUMP", fake_dump(tmp_path))
    out = tmp_path / "r.fa"
    n = pipeline.fetch("SRRX", 12, str(out), ranges=pipeline.spot_ranges(12, 100, 3), wait=0)
    ids = [l.strip() for l in open(out) if l.startswith(">")]
    assert n == 12 and ids[0] == ">r1" and ">r34" in ids and ">r67" in ids


# --- reads: stop codons, GC, length, missing paralogs ----------------------------------------

def test_stop_codons_flag_pseudogene(tmp_path):
    f = basic(tmp_path)
    rows = reads("a", "a1", 60) + reads("b", "b1", 60) + reads("c", "c1", 30, stop=True)
    s = summary(tmp_path, f, rows)
    assert s["stop_fraction_candidate"] == 1.0 and s["stop_fraction_controls"] == 0.0
    assert "possible_pseudogene" in s["flags"]


def test_no_pseudogene_flag_with_sequencing_error_background(tmp_path):
    f = basic(tmp_path)
    rows = reads("a", "a1", 60) + reads("b", "b1", 60) + reads("c", "c1", 58) + reads("x", "c1", 1, stop=True)
    assert "possible_pseudogene" not in summary(tmp_path, f, rows)["flags"]


def test_read_gc_and_short_reads(tmp_path):
    f = basic(tmp_path)
    short = "GC" * 40                                    # 80 nt, GC 1.0
    rows = reads("a", "a1", 60, seq=short) + reads("b", "b1", 60, seq=short) + reads("c", "c1", 60, seq=short)
    s = summary(tmp_path, f, rows)
    assert s["mean_read_length"] == 80 and s["gc_controls"] == 1.0
    assert "short_reads" in s["flags"]


def test_old_tables_still_work(tmp_path):
    f = basic(tmp_path)
    rows = [r[:7] for r in reads("a", "a1", 60) + reads("b", "b1", 60) + reads("c", "c1", 60)]
    s = summary(tmp_path, f, rows)
    assert s["status"] == stats.PRESENT and math.isnan(s["stop_fraction_candidate"])


def test_discarded_paralog_fragments_alone_are_not_flagged(tmp_path):
    # a species without the gene: fragments of an unlisted paralog are discarded and the absence is safe
    f = basic(tmp_path)
    rows = reads("a", "a1", 60) + reads("b", "b1", 60)
    rows += reads("d", "c1", 20, pid=62, start=10, step=1, span=20)        # paralog fragments in one domain
    s = summary(tmp_path, f, rows)
    assert s["status"] == stats.ABSENT and s["identity_only_discarded"] == 20
    assert "possible_missing_paralog" not in s["flags"]


def test_missing_paralog_when_kept_reads_are_bunched_too(tmp_path):
    f = basic(tmp_path)
    rows = reads("a", "a1", 60) + reads("b", "b1", 60)
    rows += reads("d", "c1", 20, pid=62, start=10, step=1, span=20)        # discarded, low identity
    rows += reads("k", "c1", 20, pid=80, start=10, step=1, span=20)        # kept, same domain
    s = summary(tmp_path, f, rows)
    assert "possible_missing_paralog" in s["flags"]


# --- statistics ------------------------------------------------------------------------------

def test_conditional_limits_are_wider_than_fixed():
    lo_f, hi_f = [v / 35.6 for v in stats.garwood(28)]
    lo_c, hi_c, p = stats.conditional(28, 159, 357, 1596)
    assert lo_c < lo_f and hi_c > hi_f and 0 < p < 1


def test_conditional_exact_coverage():
    # with equal lengths and true ratio 1, the 95% limits cover 1 in about 95% of draws
    import random
    rng = random.Random(3)
    hit = 0
    for _ in range(400):
        x = sum(rng.random() < 0.5 for _ in range(40))
        lo, hi, _ = stats.conditional(x, 40 - x, 100, 100)
        hit += lo <= 1 <= hi
    assert hit / 400 > 0.93


def test_copy_number_states():
    assert stats.copy_number_states(0, 300, 100, 200)[0] == 0.0
    assert stats.copy_number_states(75, 300, 100, 200)[0] == 0.5
    assert stats.copy_number_states(150, 300, 100, 200)[0] == 1.0
    assert stats.copy_number_states(300, 300, 100, 200)[0] == 2.0
    best, lr, _ = stats.copy_number_states(1, 2, 100, 200)
    assert lr < 1                                           # few reads: weak support


def test_duplicated_controls_declared(tmp_path):
    f = basic(tmp_path, ctrl_copies=2)                     # controls kept in two copies after a duplication
    rows = reads("a", "a1", 120) + reads("b", "b1", 120) + reads("c", "c1", 60)
    s = summary(tmp_path, f, rows)
    assert s["ratio"] == pytest.approx(1.0) and s["copy_state"] == 1.0


# --- sex-linked genes --------------------------------------------------------------------------

def test_z_linked_candidate_in_female(tmp_path):
    f = basic(tmp_path, cand_loc="Z")
    rows = reads("a", "a1", 100) + reads("b", "b1", 100) + reads("c", "c1", 50)     # half dosage
    assert summary(tmp_path, f, rows, sex="female")["ratio"] == pytest.approx(1.0)
    s = summary(tmp_path, f, rows)
    assert s["ratio"] == pytest.approx(0.5) and "sex_unknown_for_sex_linked_gene" in s["flags"]


def test_w_gene_in_male_is_uninformative(tmp_path):
    f = basic(tmp_path, cand_loc="W")
    rows = reads("a", "a1", 100) + reads("b", "b1", 100)
    s = summary(tmp_path, f, rows, sex="male")
    assert s["status"] == stats.UNINFORMATIVE and "not_expected_in_this_sex" in s["flags"]


def test_sex_from_markers(tmp_path):
    g, p = fam_files(tmp_path, [("cand", "candidate", 100, "", "W"), ("a", "control", 100, "", ""),
                                ("b", "control", 100, "", ""), ("z", "marker", 100, "", "Z")],
                     [("c1", "cand"), ("a1", "a"), ("b1", "b"), ("z1", "z")])
    f = family.load(g, p)
    female = reads("a", "a1", 100) + reads("b", "b1", 100) + reads("z", "z1", 50) + reads("c", "c1", 50)
    s = summary(tmp_path, f, female)
    assert s["sex"] == "female" and s["sex_source"] == "markers" and s["ratio"] == pytest.approx(1.0)
    male = reads("a", "a1", 100) + reads("b", "b1", 100) + reads("z", "z1", 100)
    s = summary(tmp_path, f, male, )
    assert s["sex"] == "male" and "not_expected_in_this_sex" in s["flags"]


def test_marker_needs_sex_location(tmp_path):
    g, p = fam_files(tmp_path, [("cand", "candidate", 100, "", ""), ("a", "control", 100, "", ""),
                                ("m", "marker", 100, "", "")], [("c1", "cand"), ("a1", "a"), ("m1", "m")])
    with pytest.raises(ValueError):
        family.load(g, p)


# --- GC and panel --------------------------------------------------------------------------

def gc_seq(gc, n=150):
    k = round(gc * n)
    return "G" * k + "A" * (n - k)


def panel_family(tmp_path, n=40):
    f = basic(tmp_path)
    t = tmp_path / "panel.tsv"
    t.write_text("protein_id\tgene\tlength_aa\n" + "".join(f"pp{i}\tpanel:{i}\t200\n" for i in range(n)))
    return family.load_panel(t, f)


def panel_reads(n=40, rate_of_gc=lambda gc: 1.0, per=60):
    rows = []
    for i in range(n):
        gc = 0.30 + 0.40 * i / (n - 1)
        rows += reads(f"P{i}_", f"pp{i}", round(per * 2 * rate_of_gc(gc)), seq=gc_seq(gc))
    return rows


def test_panel_agrees_with_single_copy_controls(tmp_path):
    f = panel_family(tmp_path)
    rows = panel_reads() + reads("a", "a1", 60) + reads("b", "b1", 60) + reads("c", "c1", 60)
    s = summary(tmp_path, f, rows)
    assert s["control_vs_panel"] == pytest.approx(1.0, abs=0.05)
    assert s["expected_panel"] == pytest.approx(60, abs=3)
    assert "controls_copy_number" not in s["flags"] and "panel_disagrees" not in s["flags"]


def test_panel_reveals_duplicated_controls(tmp_path):
    f = panel_family(tmp_path)
    rows = panel_reads() + reads("a", "a1", 120) + reads("b", "b1", 120) + reads("c", "c1", 60)
    s = summary(tmp_path, f, rows)
    assert s["control_vs_panel"] == pytest.approx(2.0, abs=0.1) and "controls_copy_number" in s["flags"]


def test_gc_bias_from_panel_curve(tmp_path):
    f = panel_family(tmp_path)
    falloff = lambda gc: 1.0 if gc < 0.55 else 0.4      # GC-rich DNA under-covered
    rows = panel_reads(rate_of_gc=falloff) + reads("a", "a1", 60, seq=gc_seq(0.42)) + \
        reads("b", "b1", 60, seq=gc_seq(0.42)) + reads("c", "c1", 24, seq=gc_seq(0.68))
    cds = {"cand": 0.68, "a": 0.42, "b": 0.42}
    s = summary(tmp_path, f, rows, cds_gc=cds)
    assert "gc_mismatch" in s["flags"] and "gc_bias" in s["flags"]
    assert s["gc_factor"] == pytest.approx(0.4, abs=0.1)
    assert s["ratio"] == pytest.approx(0.4, abs=0.05)
    corrected = summary(tmp_path, f, rows, cds_gc=cds, gc_correct=True)
    assert corrected["ratio"] == pytest.approx(1.0, abs=0.2)


def test_curve_at_outside_range():
    c = [(0.3, 1.0), (0.5, 0.8)]
    assert extras.curve_at(c, 0.4) == pytest.approx(0.9)
    assert math.isnan(extras.curve_at(c, 0.7))


def test_panel_command(tmp_path):
    faa = tmp_path / "busco.faa"
    faa.write_text("".join(f">g{i}_1\n{'M' * (250 + i)}\n>g{i}_2\n{'M' * 400}\n" for i in range(10)) +
                   ">tiny_1\nMMM\n")
    out = tmp_path / "panel.faa"
    cli.main(["panel", "--faa", str(faa), "--n", "5", "--variant-sep", "_", "--out", str(out)])
    ids = [l[1:].strip() for l in open(out) if l.startswith(">")]
    assert len(ids) == 5 and all(i.startswith("g") for i in ids)


# --- selftest: divergence per gene and across sites ------------------------------------------

def test_gene_divergence_and_gamma(tmp_path):
    f = basic(tmp_path)
    seqs = {"c1": "A" * 2000, "a1": "A" * 2000, "b1": "A" * 2000}
    fr = selftest.fragments(f, seqs, length=2000, step=2000, gene_divergence={"cand": 0.3}, divergence=0.05)
    changed = {i.split("::")[0]: sum(x != "A" for x in s) / 2000 for i, s in fr}
    assert changed["cand"] > 0.25 and changed["a"] < 0.08
    fr_g = selftest.fragments(f, {"c1": "A" * 4000, "a1": "A", "b1": "A"}, length=4000, step=4000,
                              divergence=0.2, alpha=0.5, seed=4)
    s_ = fr_g[0][1]
    windows = [sum(x != "A" for x in s_[i:i + 200]) for i in range(0, 4000, 200)]
    assert max(windows) - min(windows) > 30                    # fast and slow regions


# --- DIAMOND with the new output fields -----------------------------------------------------

@pytest.mark.skipif(not HAS_DIAMOND, reason="DIAMOND not installed")
def test_diamond_fields(tmp_path):
    prot = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAE".replace("U", "S")
    (tmp_path / "f.faa").write_text(f">P1\n{prot}\n")
    codon = {"M": "ATG", "K": "AAA", "T": "ACC", "A": "GCT", "Y": "TAT", "I": "ATT", "Q": "CAA", "R": "CGT",
             "S": "TCT", "F": "TTT", "V": "GTT", "H": "CAT", "L": "CTG", "E": "GAA", "G": "GGT", "P": "CCT",
             "D": "GAT", "N": "AAT", "W": "TGG", "C": "TGT"}
    dna = "".join(codon[c] for c in prot[:60])
    read = dna[:90] + "TAA" + dna[93:180]
    (tmp_path / "r.fa").write_text(f">r1\n{dna[:180]}\n>r2\n{read}\n")
    db = pipeline.makedb(str(tmp_path / "f.faa"), str(tmp_path / "f"))
    out = pipeline.search(str(tmp_path / "r.fa"), db, str(tmp_path / "h.tsv"), threads=1)
    rows = [l.rstrip("\n").split("\t") for l in open(out)]
    assert len(rows) == 2 and all(len(r) == 9 for r in rows)
    by = {r[0]: r for r in rows}
    assert "*" not in by["r1"][7] and "*" in by["r2"][7] and len(by["r1"][8]) == 180


@pytest.mark.skipif(not HAS_DIAMOND, reason="DIAMOND not installed")
def test_makedb_with_panel(tmp_path):
    (tmp_path / "f.faa").write_text(">F1\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEV\n")
    (tmp_path / "p.faa").write_text(">B1\nMSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVG\n>B2\nMADEEKLPPGWEKRMSRSSGRVYYFNHITNA\n")
    pipeline.makedb(str(tmp_path / "f.faa"), str(tmp_path / "fam"), panel_faa=str(tmp_path / "p.faa"))
    t = (tmp_path / "fam.panel.tsv").read_text().splitlines()
    assert t[0] == "protein_id\tgene\tlength_aa" and t[1].startswith("B1\tpanel:B1\t")
    assert (tmp_path / "fam.dmnd").exists()
