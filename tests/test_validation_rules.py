"""Rules added after the X-Y gametolog validation: assignable threshold, leakage, aligned runs."""
from readloss import assign, cli, diagnostics, family, pipeline, stats


def fam(tmp_path, assignable=None, leak_cand=0):
    (tmp_path / "g.tsv").write_text("gene\trole\tlength_aa\ncand\tcandidate\t100\nctrl\tcontrol\t100\n")
    (tmp_path / "p.tsv").write_text("protein_id\tgene\npC\tcand\npK\tctrl\n")
    f = family.load(tmp_path / "g.tsv", tmp_path / "p.tsv")
    if assignable is not None:
        (tmp_path / "st.tsv").write_text(
            "gene\trole\tfragments\tassigned_candidate\tassignable_fraction\n"
            f"cand\tcandidate\t100\t90\t{assignable}\nctrl\tcontrol\t100\t{leak_cand}\t0.9\n")
        family.load_assignable(tmp_path / "st.tsv", f)
    return f


def counts(cand, ctrl):
    c = assign.Counts()
    c.by_gene.update({"cand": cand, "ctrl": ctrl})
    return c


def test_threshold_is_one_half():
    assert diagnostics.ASSIGNABLE_MIN == 0.5


def test_candidate_below_one_half_is_uninformative(tmp_path):
    f = fam(tmp_path, assignable=0.37)
    s = assign.summarise(counts(20, 50), f, rule=assign.Rule())
    assert s["status"] == stats.UNINFORMATIVE and "low_assignable_fraction" in s["flags"]


def test_leak_is_read_from_selftest(tmp_path):
    f = fam(tmp_path, assignable=0.9, leak_cand=10)
    assert abs(f.leak["ctrl"] - 0.10) < 1e-9


def test_reads_explained_by_leakage_are_not_a_presence(tmp_path):
    f = fam(tmp_path, assignable=0.9, leak_cand=10)        # 10% of control fragments go to the candidate
    s = assign.summarise(counts(5, 45), f, rule=assign.Rule())   # 45 / 0.9 x 0.1 = 5 leaked reads expected
    assert abs(float(s["expected_if_absent"]) - 5.0) < 1e-6
    assert s["status"] == stats.UNINFORMATIVE and "leakage_explains_candidate" in s["flags"]


def test_reads_well_above_leakage_stay_present(tmp_path):
    f = fam(tmp_path, assignable=0.9, leak_cand=1)
    s = assign.summarise(counts(45, 45), f, rule=assign.Rule())
    assert s["status"] == stats.PRESENT and "leakage_explains_candidate" not in s["flags"]
    assert float(s["p_leakage"]) < 1e-6


def test_no_leak_columns_keep_old_behaviour(tmp_path):
    f = fam(tmp_path)
    s = assign.summarise(counts(45, 45), f, rule=assign.Rule())
    assert s["expected_if_absent"] == "" and s["status"] == stats.PRESENT


ALIGNED = """<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE><RUN_SET><RUN>
<SRAFiles><SRAFile filename="x.final.cram" semantic_name="cram" supertype="Original"/></SRAFiles>
<AlignInfo path="p" cnt="2"><Alignment name="chr1" seqid="NC_1"/></AlignInfo>
</RUN></RUN_SET></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>"""
FASTQ = """<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE><RUN_SET><RUN>
<SRAFiles><SRAFile filename="a.R1.fq.gz" semantic_name="fastq" supertype="Original"/></SRAFiles>
</RUN></RUN_SET></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>"""


def test_aligned_submission_is_detected():
    assert pipeline.parse_aligned(ALIGNED)
    assert not pipeline.parse_aligned(FASTQ)
    assert not pipeline.parse_aligned("not xml")
    bam_only = FASTQ.replace('a.R1.fq.gz" semantic_name="fastq"', 'x.bam" semantic_name="bam"')
    assert pipeline.parse_aligned(bam_only)


def test_only_the_requested_run_is_read():
    two = ("<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE><RUN_SET>"
           '<RUN accession="ERR1"><SRAFiles><SRAFile semantic_name="cram" supertype="Original"/></SRAFiles>'
           '<AlignInfo/></RUN>'
           '<RUN accession="ERR2"><SRAFiles><SRAFile semantic_name="fastq" supertype="Original"/></SRAFiles></RUN>'
           "</RUN_SET></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>")
    assert pipeline.parse_aligned(two, "ERR1") and not pipeline.parse_aligned(two, "ERR2")
    assert not pipeline.parse_aligned(two, "ERR3")


def test_part_of_a_sorted_run_is_refused():
    info = {"Aligned": "yes", "spots": "100000000"}
    assert cli._sorted_sample(info, 30_000_000)
    assert not cli._sorted_sample(info, 100_000_000)          # the whole run is a fair sample
    assert not cli._sorted_sample({"Aligned": "no", "spots": "100000000"}, 30_000_000)
    assert not cli._sorted_sample({}, 30_000_000)


def test_control_identity_ignores_distant_family_reads(tmp_path):
    f = fam(tmp_path)
    c = counts(10, 40)
    c.candidate_identity = [86.0] * 10
    c.identity["ctrl"] = [99.0] * 20 + [40.0] * 20               # half of them from distant family members
    s = assign.summarise(c, f, rule=assign.Rule())
    assert abs(s["identity_gap"] - 13.0) < 1e-6 and "divergent_candidate" in s["flags"]


def test_presence_from_an_unlisted_paralog_is_uninformative(tmp_path):
    f = fam(tmp_path)
    c = counts(40, 40)
    c.candidate_identity = [80.0] * 40                           # less similar than the controls ...
    c.identity["ctrl"] = [98.0] * 40
    c.candidate_spans = [(10, 40)] * 40                          # ... and bunched in one region
    c.identity_only_spans = [(10, 40)] * 6
    s = assign.summarise(c, f, rule=assign.Rule())
    assert "possible_missing_paralog" in s["flags"] and "divergent_candidate" in s["flags"]
    assert s["status"] == stats.UNINFORMATIVE and "presence_from_unlisted_paralog" in s["flags"]


def test_nucleotide_reads_of_another_species_make_a_low_ratio_uninformative(tmp_path):
    f = fam(tmp_path)
    nt = assign.Rule(min_aln=30, margin=5, min_identity=95.0)
    c = counts(5, 45)
    c.candidate_identity = [96.5] * 5                            # a horse read against a cattle sequence
    c.identity["ctrl"] = [99.8] * 45
    s = assign.summarise(c, f, rule=nt)
    assert "candidate_reference_divergent" in s["flags"] and s["status"] == stats.UNINFORMATIVE
    c.candidate_identity = [99.7] * 5                            # own-species sequence: the low ratio stands
    s = assign.summarise(c, f, rule=nt)
    assert "candidate_reference_divergent" not in s["flags"] and s["status"] != stats.UNINFORMATIVE


def test_divergent_reference_does_not_apply_to_protein_mode(tmp_path):
    f = fam(tmp_path)
    c = counts(5, 45)
    c.candidate_identity = [90.0] * 5
    c.identity["ctrl"] = [92.0] * 45
    s = assign.summarise(c, f, rule=assign.Rule())
    assert "candidate_reference_divergent" not in s["flags"]


def test_selftest_rerun_replaces_its_table(tmp_path):
    import os
    ex = os.path.join(os.path.dirname(__file__), "..", "ejemplos", "dcx")
    sim = tmp_path / "sim.tsv"
    sim.write_text("DCX::1\tCaraur_LOC113084161_113084161\t99\t50\t100\n"
                   "dclk1a::1\tCaraur_LOC113084161_113084161\t60\t50\t40\n")
    out = tmp_path / "st.tsv"
    args = ["selftest", "--hits", str(sim), "--genes", os.path.join(ex, "genes.tsv"),
            "--proteins", os.path.join(ex, "proteins.tsv"), "--out", str(out)]
    cli.main(args)
    n = len(out.read_text().splitlines())
    cli.main(args)
    assert len(out.read_text().splitlines()) == n
