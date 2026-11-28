"""Command line interface.

readloss makedb     build the DIAMOND database of the gene family
readloss test       test one run (SRA accession, local FASTA or existing hit table)
readloss batch      test many runs, resumable, with automatic deeper repeat of intermediates
readloss subsample  ratio and limits against read depth, from one hit table
readloss calibrate  sweep the assignment rule on samples of known status
readloss sensitive  slower tBLASTn + BLASTx test for divergent queries or few reads
readloss depth      reads needed, from a short pilot run or an existing hit table
readloss selftest   fraction of reads of each gene the rule can assign (run once per family)
readloss panel      choose a panel of single-copy genes from a protein FASTA (for makedb --panel)
readloss report     summary of a results table: calls, warnings and an optional figure
readloss plot       figures to decide what to trust: selftest, power, coverage, margins, identity
"""
import argparse
import csv
import os
import re
import sys
import tempfile

from . import __version__, assign, family as fam, pipeline, selftest, sensitive, stats

FMT = {"expected": "{:.1f}", "ratio": "{:.3f}", "ci_low": "{:.3f}", "ci_high": "{:.3f}",
       "p_le_obs_if_present": "{:.2g}", "mean_identity_candidate": "{:.1f}", "control_fold": "{:.2f}",
       "control_p": "{:.2g}", "ratio_loo_min": "{:.3f}", "ratio_loo_max": "{:.3f}", "p_miss_one_copy": "{:.2g}",
       "discarded_fraction": "{:.2f}", "discarded_identity": "{:.1f}", "identity_gap": "{:.1f}", "breadth": "{:.2f}",
       "breadth_expected": "{:.2f}", "largest_gap_aa": "{:.0f}", "assignable_candidate": "{:.3f}",
       "assignable_fraction": "{:.3f}", "ci_low_cond": "{:.3f}", "ci_high_cond": "{:.3f}",
       "p_cond_if_present": "{:.2g}", "copy_state_log10_lr": "{:.1f}", "stop_fraction_candidate": "{:.3f}",
       "stop_fraction_controls": "{:.3f}", "gc_candidate": "{:.3f}", "gc_controls": "{:.3f}",
       "gc_factor": "{:.2f}", "expected_panel": "{:.1f}", "control_vs_panel": "{:.2f}",
       "mean_read_length": "{:.0f}", "marker_ratio": "{:.2f}", "mate_length_metadata": "{:.0f}",
       "expected_if_absent": "{:.1f}", "p_leakage": "{:.2g}"}


def _fmt(row):
    return {k: (FMT[k].format(v) if k in FMT and isinstance(v, float) else v) for k, v in row.items()}


def _write(rows, path):
    """Write rows as TSV to stdout or append to ``path`` (header only for a new file).

    When appending, the columns of the existing header are kept; a column that the new rows add
    is a sign of a bug, so it is reported instead of being dropped silently.
    """
    rows = [_fmt(r) for r in rows]
    if not rows:
        return
    new = not path or not os.path.exists(path) or os.path.getsize(path) == 0
    names = list(rows[0])
    for r in rows[1:]:
        names += [k for k in r if k not in names]
    if not new:
        with open(path) as fh:
            head = fh.readline().rstrip("\n").split("\t")
        extra = [k for k in names if k not in head]
        if extra:
            print(f"warning: columns not in {path} were left out: {', '.join(extra)}", file=sys.stderr)
        names = head
    out = open(path, "a", newline="") if path else sys.stdout
    w = csv.DictWriter(out, fieldnames=names, delimiter="\t", extrasaction="ignore")
    if new:
        w.writeheader()
    w.writerows(rows)
    if path:
        out.close()


def _common(p):
    p.add_argument("--genes", required=True, help="genes.tsv (gene, role, length_aa)")
    p.add_argument("--proteins", required=True, help="proteins.tsv (protein_id, gene)")
    p.add_argument("--min-aln", type=int, default=None,
                   help="minimum alignment length in aa [30, or 60%% of the read length when the SRA reports "
                        "reads shorter than 100 nt]")
    p.add_argument("--copies", help="copies per haploid genome for this specimen, overriding genes.tsv, as "
                                    "gene=n,gene=n (a runs.tsv column copies does the same per run)")
    p.add_argument("--margin", type=float, default=5.0, help="bits over the best hit to any other gene [5]")
    p.add_argument("--min-identity", type=float, default=None,
                   help="identity of counted candidate reads, percent [75; 95 in nucleotide mode]")
    p.add_argument("--nt", action="store_true",
                   help="nucleotide mode: --db is a FASTA of coding sequences indexed by makedb --nt, reads are "
                        "aligned with bwa-mem2 -a and --margin is in alignment-score points (5 = one mismatch); "
                        "for homeologs and recent duplicates that encode the same protein")
    p.add_argument("--rule", choices=[assign.SYMMETRIC, assign.CENSUS], default=assign.SYMMETRIC,
                   help="symmetric: candidate and controls pass the same margin and identity, ties discarded; "
                        "census: rule of the doublecortin census (version 0.1.0), margin only for the "
                        "candidate [symmetric]")
    p.add_argument("--assignable", help="table from readloss selftest; corrects the expected count for reads "
                                        "the rule cannot assign")
    p.add_argument("--sex", help="sex of the specimen (female, male, F, M ...), for genes with a sex-linked "
                                 "location in genes.tsv [from the SRA metadata or the markers]")
    p.add_argument("--cds", help="FASTA of the coding sequences of candidate and controls (ids = gene or "
                                 "protein ids), for the GC diagnostics")
    p.add_argument("--panel", help="panel.tsv written by makedb --panel; adds the panel expectation and the "
                                   "GC curve")
    p.add_argument("--family-fasta", help="the family FASTA; read positions are mapped to one reference protein "
                                          "per gene, so coverage and breadth are exact")
    p.add_argument("--gc-correct", action="store_true",
                   help="correct the expected count with the GC curve of the panel (needs --panel and --cds)")
    p.add_argument("--absent-max", type=float, default=0.10, help="ratio scored as absence [0.10]")
    p.add_argument("--present-min", type=float, default=0.30, help="ratio scored as presence [0.30]")
    p.add_argument("--min-expected", type=float, default=10.0, help="expected reads below this are uninformative [10]")


def _search_opts(p):
    p.add_argument("--db", help="DIAMOND database of the family (from readloss makedb)")
    p.add_argument("--spots", type=int, default=15_000_000, help="spots to download first [15000000]")
    p.add_argument("--spots-max", type=int, default=40_000_000,
                   help="spots for the repeat of intermediate results; 0 disables it [40000000]")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--workdir", default=None, help="folder for reads and hit tables [temporary]")
    p.add_argument("--keep-reads", action="store_true", help="do not delete downloaded reads")
    p.add_argument("--stream", action="store_true",
                   help="pipe the download straight into DIAMOND; no reads are written to disk")
    p.add_argument("--allow-short", action="store_true", help="accept runs with fewer spots than requested")
    p.add_argument("--metadata", choices=["check", "skip"], default="check",
                   help="read the SRA metadata of the run: library type, sex, read length [check]")
    p.add_argument("--blocks", type=int, default=1,
                   help="take the spots in this many blocks spread along the run instead of from its start [1]")


def _cds_gc(path, f):
    """GC of each gene from a CDS FASTA whose ids are gene labels or protein ids of proteins.tsv."""
    if not path:
        return None
    from . import selftest as st
    from .assign import _gc
    by = {}
    for k, s_ in st.read_fasta(path).items():
        g = k if k in f.lengths else f.protein_gene.get(k)
        if g is not None:
            by.setdefault(g, []).append(_gc(s_))
    return {g: sum(v) / len(v) for g, v in by.items()}


def parse_copies(text):
    """gene=n,gene=n (also ; or spaces) -> {gene: n}."""
    out = {}
    for item in re.split(r"[,;\s]+", (text or "").strip()):
        if not item:
            continue
        g, _, n = item.partition("=")
        if not n:
            raise ValueError(f"copies {item!r}: use gene=n")
        out[g.strip()] = float(n)
        if out[g.strip()] <= 0:
            raise ValueError(f"copies of {g!r} must be positive")
    return out


def _with_copies(f, text):
    """A copy of the family with the copy numbers of one specimen (lineage with a genome duplication)."""
    import dataclasses
    over = parse_copies(text)
    unknown = [g for g in over if g not in f.lengths]
    if unknown:
        raise ValueError(f"copies given for unknown genes: {', '.join(unknown)}")
    return dataclasses.replace(f, copies={**f.copies, **over})


def _short_read_rule(rule, a, info):
    """Shorter minimum alignment for short reads, unless --min-aln was given."""
    import dataclasses
    ml = pipeline.mate_length(info) if info else float("nan")
    if a.min_aln is None and ml == ml and ml < assign_short():
        n = max(15, min(30, int(round(0.6 * ml / 3))))
        return dataclasses.replace(rule, min_aln=n)
    return rule


def _nt(a):
    return bool(getattr(a, "nt", False) or getattr(a, "sam", None))


def _min_identity(a):
    if a.min_identity is not None:
        return a.min_identity
    return 95.0 if _nt(a) else 75.0


def _search_fn(a):
    return pipeline.search_nt if _nt(a) else pipeline.search


def _setup(a):
    f = fam.load(a.genes, a.proteins)
    if getattr(a, "assignable", None):
        fam.load_assignable(a.assignable, f)
    if getattr(a, "panel", None):
        fam.load_panel(a.panel, f)
    if getattr(a, "family_fasta", None):
        fam.load_coords(a.family_fasta, f)
    rule = assign.Rule(a.min_aln if a.min_aln is not None else 30, a.margin, _min_identity(a), a.rule)
    th = stats.Thresholds(a.absent_max, a.present_min, a.min_expected)
    a._cds_gc = _cds_gc(getattr(a, "cds", None), f)
    if getattr(a, "copies", None):
        f = _with_copies(f, a.copies)
    if getattr(a, "gc_correct", False) and not (getattr(a, "panel", None) and a._cds_gc):
        sys.exit("--gc-correct needs --panel and --cds")
    return f, rule, th


def _from_hits(hits_tsv, f, rule, th, a=None, sex=None, reads=None):
    info = {}
    hits = assign.read_hits(hits_tsv, f, rule, info=info)
    sex = sex or (getattr(a, "sex", None) if a else None) or "unknown"
    return assign.summarise(assign.assign(hits, f, rule, reads=reads, info=info), f, th, rule=rule, sex=sex,
                            cds_gc=getattr(a, "_cds_gc", None) if a else None,
                            gc_correct=getattr(a, "gc_correct", False) if a else False)


def _metadata(run, a):
    """SRA metadata of a run, or {} when skipped or unreachable (noted in the log)."""
    if getattr(a, "metadata", "skip") == "skip":
        return {}
    try:
        info = pipeline.run_info(run)
    except Exception as e:
        print(f"warning: {e}", file=sys.stderr)
        return {}
    try:
        desc, aligned = pipeline.library_package(run)
        info["LibraryDescription"], info["Aligned"] = desc, "yes" if aligned else "no"
    except Exception:
        info["LibraryDescription"] = ""
    return info


def _sorted_sample(info, spots):
    """True when the run is stored sorted by position and ``spots`` is only part of it."""
    total = int(float((info or {}).get("spots") or 0))
    return (info or {}).get("Aligned") == "yes" and (not total or spots < 0.95 * total)


def _apply_metadata(row, info):
    if not info:
        return row
    issues = pipeline.library_issues(info)
    sex_meta = info.get("Sex", "")
    row.update({"library_strategy": info.get("LibraryStrategy", ""),
                "library_selection": info.get("LibrarySelection", ""),
                "library_source": info.get("LibrarySource", ""),
                "library_layout": info.get("LibraryLayout", ""),
                "mate_length_metadata": pipeline.mate_length(info),
                "sex_metadata": sex_meta, "organism": info.get("ScientificName", ""),
                "pcr_free": pipeline.pcr_free(info.get("LibraryDescription", ""))})
    flags = [x for x in (row.get("flags") or "").split(",") if x]
    row["library_issues"] = ""
    if issues:
        flags.append("library_not_wgs")
        row["status"] = stats.UNINFORMATIVE
        row["library_issues"] = "; ".join(issues)
    if "pool" in sex_meta.lower() or "pool" in (info.get("SampleName") or "").lower():
        flags.append("pooled_sample")
    ml = pipeline.mate_length(info)
    if ml == ml and ml < assign_short() and "short_reads" not in flags:
        flags.append("short_reads")
    row["flags"] = ",".join(flags)
    return row


def assign_short():
    from .diagnostics import SHORT_READ_NT
    return SHORT_READ_NT


def _one_pass(run, spots, a, f, rule, th, wd, info=None, sex=None):
    fasta = os.path.join(wd, f"{run}.{spots}.fasta")
    tsv = os.path.join(wd, f"{run}.{spots}.hits.tsv")
    spots_file = tsv + ".spots"
    if not os.path.exists(tsv):
        log = os.path.join(wd, "download.log")
        total = int((info or {}).get("spots") or 0)
        ranges = pipeline.spot_ranges(spots, total, getattr(a, "blocks", 1)) \
            if getattr(a, "blocks", 1) > 1 and total else None
        if a.stream:
            got = pipeline.fetch_search_stream(run, spots, a.db, tsv, threads=a.threads, log=log,
                                               allow_short=a.allow_short, ranges=ranges, search_fn=_search_fn(a))
        else:
            got = pipeline.fetch(run, spots, fasta, log=log, allow_short=a.allow_short, ranges=ranges)
            _search_fn(a)(fasta, a.db, tsv, threads=a.threads)
            if not a.keep_reads:
                os.remove(fasta)
        with open(spots_file, "w") as fh:
            fh.write(f"{got}\n")
    else:
        got = int(open(spots_file).read()) if os.path.exists(spots_file) else spots
    row = _from_hits(tsv, f, rule, th, a, sex=sex)
    row["spots"] = got
    row["hits_table"] = tsv
    return row


def _run_sex(a, info, given=None):
    """Sex for a run: given (option or runs.tsv) > SRA metadata > unknown (then markers)."""
    if given and fam.norm_sex(given) != "unknown":
        return fam.norm_sex(given)
    if getattr(a, "sex", None) and fam.norm_sex(a.sex) != "unknown":
        return fam.norm_sex(a.sex)
    return fam.norm_sex((info or {}).get("Sex", ""))


def _test_run(run, label, a, f, rule, th, wd, sex=None, copies=None, spots=None):
    info = _metadata(run, a)
    sx = _run_sex(a, info, sex)
    if copies:
        f = _with_copies(f, copies)
    rule = _short_read_rule(rule, a, info)
    n1 = int(float(spots)) if spots else a.spots
    sorted_part = _sorted_sample(info, max(n1, a.spots_max or 0))
    if info and (pipeline.library_issues(info) or sorted_part):
        # not a random whole-genome library, or a part of a position-sorted alignment: nothing to download
        first = assign.summarise(assign.Counts(), f, th, rule=rule, sex=sx)
        first.update({"spots": 0, "hits_table": ""})
        row = {"label": label, "run": run}
        row.update(first)
        row.update({"repeated_deeper": False, "first_spots": "", "first_observed": "", "first_expected": "",
                    "first_ratio": "", "min_aln_used": rule.min_aln,
                    "copies_used": ",".join(f"{g}={n:g}" for g, n in sorted(f.copies.items()))})
        row = _apply_metadata(row, info)
        if sorted_part:
            row["status"] = stats.UNINFORMATIVE
            row["flags"] = ",".join([x for x in (row.get("flags") or "").split(",") if x] + ["aligned_sorted_run"])
            row["library_issues"] = "; ".join([x for x in [row.get("library_issues")] if x] +
                                              ["aligned submission sorted by position: test the whole run "
                                               "(--spots at least the run size) or the original FASTQ"])
        return row
    first = _one_pass(run, n1, a, f, rule, th, wd, info, sx)
    final, repeat = first, False
    if first["status"] == stats.INTERMEDIATE and a.spots_max and a.spots_max > n1:
        final, repeat = _one_pass(run, a.spots_max, a, f, rule, th, wd, info, sx), True
    row = {"label": label, "run": run}
    row.update(final)
    row["repeated_deeper"] = repeat
    row.update({"first_spots": first["spots"] if repeat else "", "first_observed": first["observed"] if repeat else "",
                "first_expected": first["expected"] if repeat else "", "first_ratio": first["ratio"] if repeat else ""})
    if info and row.get("sex_source") == "given" and not (sex or getattr(a, "sex", None)):
        row["sex_source"] = "metadata"
    row["min_aln_used"] = rule.min_aln
    row["copies_used"] = ",".join(f"{g}={n:g}" for g, n in sorted(f.copies.items()))
    return _apply_metadata(row, info)


def cmd_makedb(a):
    if a.nt:
        if a.panel:
            sys.exit("--panel is for the protein database; the nucleotide mode has no panel")
        print(pipeline.makedb_nt(a.family, a.out))
    else:
        print(pipeline.makedb(a.family, a.out, panel_faa=a.panel))


def cmd_panel(a):
    """Choose single-copy genes for the panel: one sequence per gene, lengths in range, random subset."""
    import random
    from . import selftest as st
    seqs = st.read_fasta(a.faa)
    groups = {}
    for k, v in seqs.items():
        g = k.rsplit(a.variant_sep, 1)[0] if a.variant_sep and a.variant_sep in k else k
        if a.min_len <= len(v) <= a.max_len and (g not in groups or len(v) > len(groups[g][1])):
            groups[g] = (k, v)
    keys = sorted(groups)
    random.Random(a.seed).shuffle(keys)
    keep = sorted(keys[:a.n])
    out = open(a.out, "w") if a.out else sys.stdout
    for g in keep:
        out.write(f">{g}\n{groups[g][1]}\n")
    if a.out:
        out.close()
    print(f"{len(keep)} genes of {len(groups)} eligible", file=sys.stderr)


def cmd_test(a):
    f, rule, th = _setup(a)
    if a.hits:
        row = {"label": a.label or os.path.basename(a.hits)}
        row.update(_from_hits(a.hits, f, rule, th, a))
    elif a.sam:
        wd = a.workdir or tempfile.mkdtemp(prefix="readloss_")
        os.makedirs(wd, exist_ok=True)
        tsv = pipeline.sam_file_to_hits(a.sam, os.path.join(wd, os.path.basename(a.sam) + ".hits.tsv"))
        row = {"label": a.label or os.path.basename(a.sam)}
        row.update(_from_hits(tsv, f, rule, th, a))
        row["hits_table"] = tsv
    elif a.reads:
        if not a.db:
            sys.exit("--db is required with --reads")
        wd = a.workdir or tempfile.mkdtemp(prefix="readloss_")
        tsv = os.path.join(wd, os.path.basename(a.reads) + ".hits.tsv")
        _search_fn(a)(a.reads, a.db, tsv, threads=a.threads)
        row = {"label": a.label or os.path.basename(a.reads)}
        row.update(_from_hits(tsv, f, rule, th, a))
        row["hits_table"] = tsv
    else:
        if not a.db:
            sys.exit("--db is required with --run")
        wd = a.workdir or tempfile.mkdtemp(prefix="readloss_")
        os.makedirs(wd, exist_ok=True)
        row = _test_run(a.run, a.label or a.run, a, f, rule, th, wd)
    _write([row], a.out)


def cmd_batch(a):
    f, rule, th = _setup(a)
    wd = a.workdir or "readloss_work"
    os.makedirs(wd, exist_ok=True)
    done = set()
    if os.path.exists(a.out):
        with open(a.out) as fh:
            done = {r["run"] for r in csv.DictReader(fh, delimiter="\t")}
    with open(a.runs) as fh:
        todo = [r for r in csv.DictReader(fh, delimiter="\t")]
    for r in todo:
        if r["run"] in done:
            continue
        try:
            row = _test_run(r["run"], r.get("label") or r["run"], a, f, rule, th, wd, sex=r.get("sex"),
                            copies=r.get("copies"), spots=(r.get("spots") or "").strip() or None)
        except Exception as e:                       # keep going, log the failure
            with open(os.path.join(wd, "failed.tsv"), "a") as fh:
                fh.write(f"{r['run']}\t{e}\n")
            continue
        _write([row], a.out)


def cmd_subsample(a):
    f, rule, th = _setup(a)
    info = {}
    hits = assign.read_hits(a.hits, f, rule, info=info)
    rows = []
    for fr in [float(x) for x in a.fractions.split(",")]:
        for rep in range(1, a.replicates + 1):
            keep = assign.subsample_ids(hits, fr, seed=rep)
            row = {"fraction": fr, "replicate": rep}
            if a.spots:
                row["spots"] = round(a.spots * fr)
            row.update(assign.summarise(assign.assign(hits, f, rule, reads=keep, info=info), f, th, rule=rule,
                                        sex=a.sex or "unknown", cds_gc=a._cds_gc, gc_correct=a.gc_correct))
            rows.append(row)
    _write(rows, a.out)


def cmd_depth(a):
    """Expected candidate reads per million spots and the depth that reaches the target."""
    f, rule, th = _setup(a)
    if a.hits:
        if not a.spots:
            sys.exit("--spots (spots behind the hit table) is required with --hits")
        row = _from_hits(a.hits, f, rule, th, a)
        spots, label = a.spots, a.label or os.path.basename(a.hits)
    else:
        if not a.db:
            sys.exit("--db is required with --run")
        a.stream = True
        info = _metadata(a.run, a)
        sorted_part = _sorted_sample(info, a.pilot_spots)
        if info and (pipeline.library_issues(info) or sorted_part):
            # not a whole-genome library, or a pilot would read one stretch of a position-sorted alignment:
            # report it without downloading anything
            out = {"label": a.label or a.run, "spots_used": 0, "expected": float("nan"),
                   "expected_per_million_spots": float("nan"), "target_expected": th.min_expected,
                   "safety": a.safety, "recommended_spots": "nan", "flags": ""}
            out = _apply_metadata(out, info)
            if sorted_part:
                out["flags"] = ",".join([x for x in out["flags"].split(",") if x] + ["aligned_sorted_run"])
            out.pop("status", None)
            _write([out], a.out)
            return
        wd = a.workdir or tempfile.mkdtemp(prefix="readloss_")
        os.makedirs(wd, exist_ok=True)
        row = _one_pass(a.run, a.pilot_spots, a, f, _short_read_rule(rule, a, info), th, wd, info,
                        _run_sex(a, info))
        spots, label = row["spots"], a.label or a.run
    rate = row["expected"] / spots * 1e6
    need = a.safety * th.min_expected / rate * 1e6 if rate > 0 else float("nan")
    out = {"label": label, "spots_used": spots, "expected": row["expected"],
           "expected_per_million_spots": round(rate, 3), "target_expected": th.min_expected,
           "safety": a.safety, "recommended_spots": int(-(-need // 1e6) * 1e6) if need == need else "nan",
           "flags": row.get("flags", "")}
    if not a.hits:
        out = _apply_metadata(out, info)
        out.pop("status", None)
        if out.get("library_issues"):
            out["recommended_spots"] = "nan"     # not a whole-genome library: no depth makes it usable
    _write([out], a.out)


FLAG_HELP = {
    "controls_disagree": "controls collect reads at different rates; one may not be single copy",
    "call_depends_on_one_control": "leaving out one control changes the status",
    "low_power_for_one_copy": "an absence call that would also be made for a gene present in one copy",
    "rule_may_discard_true_reads": "many candidate reads discarded with high identity; calibrate the rule",
    "divergent_candidate": "candidate reads much less similar than control reads; consider the sensitive mode",
    "possible_partial_loss": "reads cover less of the candidate than expected; part of the gene may be missing",
    "low_assignable_fraction": "the rule can assign few candidate reads (selftest); the status is uninformative",
    "possible_pseudogene": "many candidate reads carry stop codons; the gene may be a pseudogene or retrocopy",
    "possible_missing_paralog": "kept candidate reads bunched in one domain look like a family member missing from the database",
    "gc_mismatch": "candidate and controls differ in GC; Illumina coverage may differ between them",
    "gc_bias": "the GC curve of the panel predicts different coverage for candidate and controls",
    "controls_copy_number": "controls collect reads at a rate unlike the panel; they may not be single copy",
    "panel_disagrees": "the panel expectation gives another status than the family controls",
    "short_reads": "reads are short for 30 aa alignments; detection is weaker",
    "sex_unknown_for_sex_linked_gene": "a sex-linked gene was tested without knowing the sex of the specimen",
    "not_expected_in_this_sex": "the candidate is not expected in this sex (Y or W gene)",
    "library_not_wgs": "the SRA metadata do not describe a random whole-genome library; status set to uninformative",
    "pooled_sample": "the run pools several individuals; the ratio reflects a mixture",
    "possible_contamination": "a few candidate reads much closer to the queries than the controls: another species",
    "leakage_explains_candidate": "the candidate reads are no more than the control reads the rule misassigns "
                                  "(selftest); status set to uninformative",
    "presence_from_unlisted_paralog": "the candidate reads are bunched and less similar than the controls, like "
                                      "an unlisted family member; status set to uninformative",
    "candidate_reference_divergent": "nucleotide mode: candidate reads match their coding sequence below 98.5%, so "
                                     "it comes from another species and reads are lost; a low ratio is uninformative",
    "aligned_sorted_run": "the run was submitted as aligned reads sorted by position (BAM/CRAM); a part of it is "
                          "not a random sample of the genome, so it was not tested",
}


def cmd_report(a):
    import collections
    with open(a.results) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")]
    st = collections.Counter(r.get("status", "") for r in rows)
    fl = collections.Counter(f for r in rows for f in (r.get("flags") or "").split(",") if f)
    lines = [f"runs\t{len(rows)}"] + [f"status_{k}\t{v}" for k, v in sorted(st.items())]
    lines += [f"flag_{k}\t{v}\t{FLAG_HELP.get(k, '')}" for k, v in fl.most_common()]
    lines += ["", "label\tstatus\tratio\tflags"]
    lines += [f"{r.get('label', '')}\t{r.get('status', '')}\t{r.get('ratio', '')}\t{r.get('flags', '')}"
              for r in rows if r.get("flags") or r.get("status") in ("intermediate", "uninformative")]
    text = "\n".join(lines) + "\n"
    if a.out:
        open(a.out, "w").write(text)
    else:
        sys.stdout.write(text)
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        th = stats.Thresholds(a.absent_max, a.present_min, a.min_expected)

        def num(x):
            try:
                return float(x)
            except (TypeError, ValueError):
                return float("nan")
        pts = []
        for r in rows:
            v = num(r.get("ratio"))
            if v == v:
                lo = num(r.get("ci_low_cond") or r.get("ci_low"))
                hi = num(r.get("ci_high_cond") or r.get("ci_high"))
                pts.append((v, lo, hi, r.get("status", ""), r.get("label") or r.get("run", "")))
        pts.sort()
        left_out = len(rows) - len(pts)
        col = {"present": "#0072B2", "absent": "#D55E00", "intermediate": "#E69F00", "uninformative": "#BDBDBD"}
        n = len(pts)
        fig, ax = plt.subplots(figsize=(max(4, min(12, 0.45 * n + 2)), 3.2))
        ax.axhspan(th.absent_max, th.present_min, color="#EBEBEB", lw=0)
        ax.axhline(1, color="#9E9E9E", lw=0.6, ls=":")
        seen = set()
        for i, (v, lo, hi, st, lab) in enumerate(pts):
            c = col.get(st, "#000000")
            if lo == lo and hi == hi and hi != float("inf"):
                ax.plot([i, i], [lo, hi], color=c, lw=1, alpha=0.6)
            ax.plot(i, v, "D", ms=5 if n <= 40 else 3, color=c, mec="#000000" if st == "intermediate" else c,
                    mew=0.5, label=st if st not in seen else None)
            seen.add(st)
        if n <= 40:
            ax.set_xticks(range(n))
            ax.set_xticklabels([p_[4] for p_ in pts], rotation=60, ha="right", fontsize=7)
        else:
            ax.set_xlabel("Runs ordered by ratio")
        ax.set_xlim(-0.6, max(n - 0.4, 0.6))
        ax.set_ylabel("Observed / expected")
        if seen:
            ax.legend(frameon=False, fontsize=7, loc="upper left")
        if left_out:
            ax.set_title(f"{left_out} run(s) without a ratio not shown (uninformative)", fontsize=8, loc="right")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        fig.savefig(a.plot, bbox_inches="tight", dpi=200)


def cmd_plot(a):
    from . import plots
    out = a.out
    if a.kind == "selftest":
        if not a.selftest:
            sys.exit("plot selftest needs --selftest (the table written by readloss selftest)")
        with open(a.selftest) as fh:
            data = plots.selftest(list(csv.DictReader(fh, delimiter="\t")), out)
    elif a.kind == "power":
        if not a.results:
            sys.exit("plot power needs --results")
        with open(a.results) as fh:
            data = plots.power(list(csv.DictReader(fh, delimiter="\t")), out, a.absent_max, a.present_min,
                               a.min_expected)
    else:
        if not (a.genes and a.proteins):
            sys.exit(f"plot {a.kind} needs --genes and --proteins")
        row, hits = {}, a.hits
        if not hits:
            if not (a.results and a.run):
                sys.exit(f"plot {a.kind} needs --hits, or --results with --run (run accession or label)")
            with open(a.results) as fh:
                rows = [r for r in csv.DictReader(fh, delimiter="\t") if a.run in (r.get("run"), r.get("label"))]
            if not rows:
                sys.exit(f"run {a.run!r} is not in {a.results}")
            row = rows[0]
            hits = row.get("hits_table")
            if not hits:
                sys.exit(f"run {a.run!r} has no hit table (status {row.get('status')}, {row.get('flags')})")
        f = fam.load(a.genes, a.proteins)
        if a.assignable:
            fam.load_assignable(a.assignable, f)
        if a.family_fasta:
            fam.load_coords(a.family_fasta, f)
        elif a.kind == "coverage":
            print("warning: without --family-fasta positions are on whichever ortholog each read hit; "
                  "coverage is approximate", file=sys.stderr)
        mode = row.get("rule") or a.rule
        min_aln = int(float(row.get("min_aln_used") or 0)) or (a.min_aln if a.min_aln is not None else 30)
        rule = assign.Rule(min_aln, a.margin, _min_identity(a), mode)
        info = {}
        h = assign.read_hits(hits, f, rule, info=info)
        counts = assign.assign(h, f, rule, info=info)
        th = stats.Thresholds(a.absent_max, a.present_min, a.min_expected)
        res = assign.summarise(counts, f, th, rule=rule)
        label = row.get("label") or a.label or os.path.basename(hits)
        if a.kind == "coverage":
            data = plots.coverage(counts, f, res, out, label)
        elif a.kind == "margins":
            data = plots.margins(counts, f, rule, out, label)
        else:
            data = plots.identity(counts, f, out, label)
    print(f"{out}\n{data}")


def cmd_calibrate(a):
    f, _, th = _setup(a)
    with open(a.samples) as fh:
        samples = [r for r in csv.DictReader(fh, delimiter="\t")]
    rows = []
    for m in [float(x) for x in a.margins.split(",")]:
        for i in [float(x) for x in a.identities.split(",")]:
            rule = assign.Rule(a.min_aln if a.min_aln is not None else 30, m, i, a.rule)
            ratios = {"present": [], "absent": []}
            for s in samples:
                res = _from_hits(s["hits"], f, rule, th, a, sex=s.get("sex"))
                ratios.setdefault(s["known_status"], []).append(res["ratio"])
                rows.append({"margin": m, "min_identity": i, "label": s["label"],
                             "known_status": s["known_status"], "observed": res["observed"],
                             "expected": res["expected"], "ratio": res["ratio"], "status": res["status"]})
            if ratios["present"] and ratios["absent"]:
                gap = min(ratios["present"]) - max(ratios["absent"])
                rows.append({"margin": m, "min_identity": i, "label": "SEPARATION",
                             "known_status": "", "observed": "", "expected": "",
                             "ratio": gap, "status": "separated" if gap > 0 else "overlap"})
    _write(rows, a.out)


def cmd_selftest(a):
    f = fam.load(a.genes, a.proteins)
    rule = assign.Rule(a.min_aln if a.min_aln is not None else 30, a.margin, _min_identity(a), a.rule)
    if a.hits:
        rows = selftest.run_hits(f, a.hits, rule)
    else:
        if not a.family:
            sys.exit("selftest needs --family (protein FASTA) or --hits (simulated reads named GENE::...)")
        gd = selftest.read_divergence(a.divergence_table) if a.divergence_table else None
        if a.nt:
            rows = selftest.run_nt(f, a.family, rule, length=a.length if a.length != 50 else 150,
                                   step=a.step, error=a.error, divergence=a.divergence, seed=a.seed,
                                   threads=a.threads, workdir=a.workdir, gene_divergence=gd)
        else:
            rows = selftest.run_protein(f, a.family, rule, db=a.db, length=a.length, step=a.step,
                                        divergence=a.divergence, seed=a.seed, threads=a.threads,
                                        workdir=a.workdir, gene_divergence=gd, alpha=a.gamma_alpha)
    for r in rows:
        r["rule"] = rule.mode
    if a.out and os.path.exists(a.out):
        os.remove(a.out)        # one table per family: a rerun replaces it instead of appending duplicate genes
    _write(rows, a.out)
    cand = [r for r in rows if r["role"] == "candidate"][0]
    if cand["assignable_fraction"] < 0.2:
        print(f"warning: only {cand['assignable_fraction']:.0%} of candidate reads can be assigned; "
              "tests of this family will be uninformative", file=sys.stderr)


def cmd_sensitive(a):
    f = fam.load(a.genes, a.proteins)
    if a.assignable:
        # selftest measures the DIAMOND search; for the sensitive mode it is an approximation
        fam.load_assignable(a.assignable, f)
    if a.copies:
        f = _with_copies(f, a.copies)
    th = stats.Thresholds(a.absent_max, a.present_min, a.min_expected)
    wd = a.workdir or tempfile.mkdtemp(prefix="readloss_sens_")
    os.makedirs(wd, exist_ok=True)
    bx = os.path.join(wd, "blastx.tsv")
    spots = None
    if a.blastx:
        bx = a.blastx
    else:
        db = os.path.join(wd, "reads_db")
        spots = sensitive.build_reads_db(a.reads, db, run=a.run, spots=a.spots,
                                         log=os.path.join(wd, "download.log"))
        sensitive.run_search(db, a.queries, a.family, wd, threads=a.threads)
        if not a.keep_db:
            for x in os.listdir(wd):
                if x.startswith("reads_db"):
                    os.remove(os.path.join(wd, x))
    row = {"label": a.label or a.run or os.path.basename(a.reads or a.blastx), "mode": "sensitive"}
    if spots:
        row["spots"] = spots
    if a.rule == assign.CENSUS:
        counts = sensitive.assign_sensitive(bx, f)        # rule of the published genome tests
    else:
        rule = assign.Rule(a.min_aln if a.min_aln is not None else 30, a.margin, _min_identity(a),
                           assign.SYMMETRIC)
        counts = assign.assign(assign.read_hits(bx, f, rule), f, rule)
    row["rule"] = a.rule
    row.update(assign.summarise(counts, f, th, sex=a.sex or "unknown", cds_gc=_cds_gc(a.cds, f)))
    _write([row], a.out)


def main(argv=None):
    p = argparse.ArgumentParser(prog="readloss", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version=f"readloss {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("makedb", help="build the DIAMOND database of the family")
    s.add_argument("--family", required=True, help="FASTA with every protein of the family")
    s.add_argument("--panel", help="FASTA of single-copy panel genes (readloss panel); writes <out>.panel.tsv")
    s.add_argument("--nt", action="store_true", help="--family is a FASTA of coding sequences; index it with "
                                                     "bwa-mem2 for the nucleotide mode")
    s.add_argument("--out", required=True, help="output database name")
    s.set_defaults(fn=cmd_makedb)

    s = sub.add_parser("panel", help="choose panel genes from a protein FASTA (e.g. BUSCO ancestral proteins)")
    s.add_argument("--faa", required=True, help="protein FASTA of single-copy genes of the clade")
    s.add_argument("--n", type=int, default=200, help="genes to keep [200]")
    s.add_argument("--min-len", type=int, default=300)
    s.add_argument("--max-len", type=int, default=1500)
    s.add_argument("--variant-sep", default="", help="separator of variant suffixes to merge, e.g. _ [none]")
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--out", help="output FASTA [stdout]")
    s.set_defaults(fn=cmd_panel)

    s = sub.add_parser("test", help="test one run")
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument("--run", help="SRA run accession")
    g.add_argument("--reads", help="local FASTA of reads")
    g.add_argument("--hits", help="existing DIAMOND table (qseqid sseqid pident length bitscore)")
    g.add_argument("--sam", help="existing SAM/BAM of reads aligned to the coding sequences (bwa-mem2 or "
                                 "minimap2 with all alignments, e.g. -a); nucleotide mode")
    s.add_argument("--label")
    s.add_argument("--out", help="output TSV [stdout]")
    _common(s); _search_opts(s)
    s.set_defaults(fn=cmd_test)

    s = sub.add_parser("batch", help="test many runs")
    s.add_argument("--runs", required=True, help="TSV with columns run and label, and optionally spots "
                                                 "(from readloss depth), sex and copies")
    s.add_argument("--out", required=True, help="results TSV; finished runs are skipped on restart")
    _common(s); _search_opts(s)
    s.set_defaults(fn=cmd_batch)

    s = sub.add_parser("subsample", help="ratio against depth from one hit table")
    s.add_argument("--hits", required=True)
    s.add_argument("--fractions", default="0.01,0.02,0.05,0.1,0.2,0.5,1")
    s.add_argument("--replicates", type=int, default=10)
    s.add_argument("--spots", type=int, help="spots behind the hit table, to report depth in spots")
    s.add_argument("--out")
    _common(s)
    s.set_defaults(fn=cmd_subsample)

    s = sub.add_parser("depth", help="spots needed, from a pilot run or a hit table")
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument("--run", help="SRA run accession for a streamed pilot")
    g.add_argument("--hits", help="existing DIAMOND table (give --spots)")
    s.add_argument("--spots", type=int, help="spots behind --hits")
    s.add_argument("--pilot-spots", type=int, default=2_000_000, help="spots for the pilot [2000000]")
    s.add_argument("--safety", type=float, default=1.5,
                   help="factor over the minimum of expected reads, for the spread between runs [1.5]")
    s.add_argument("--db"); s.add_argument("--threads", type=int, default=4); s.add_argument("--workdir")
    s.add_argument("--metadata", choices=["check", "skip"], default="check",
                   help="read the SRA metadata of the run: library type, sex, read length [check]")
    s.add_argument("--blocks", type=int, default=1, help="spread the pilot spots along the run [1]")
    s.add_argument("--allow-short", action="store_true"); s.add_argument("--label"); s.add_argument("--out")
    _common(s)
    s.set_defaults(fn=cmd_depth, keep_reads=False)

    s = sub.add_parser("selftest", help="fraction of reads of each gene the rule can assign")
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument("--family", help="FASTA with every protein of the family (fragments are cut from it)")
    g.add_argument("--hits", help="hit table of reads simulated elsewhere, identifiers GENE::anything")
    s.add_argument("--db", help="DIAMOND database of the family [built from --family]")
    s.add_argument("--length", type=int, default=50,
                   help="fragment length in aa, read length / 3 [50]; in nucleotide mode, read length in nt [150]")
    s.add_argument("--error", type=float, default=0.002,
                   help="sequencing error per base of simulated reads, nucleotide mode [0.002]")
    s.add_argument("--step", type=int, default=5, help="distance between fragment starts in aa [5]")
    s.add_argument("--divergence", type=float, default=0.0,
                   help="fraction of residues substituted at random, to mimic a divergent species [0]")
    s.add_argument("--divergence-table", help="TSV gene, divergence: one rate per gene, for a candidate that "
                                              "evolves faster than its controls")
    s.add_argument("--gamma-alpha", type=float, default=0.0,
                   help="shape of gamma rates across sites (0.5-1 is typical; 0 = equal rates) [0]")
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--threads", type=int, default=4)
    s.add_argument("--workdir")
    s.add_argument("--out", help="output TSV, to pass to --assignable [stdout]")
    _common(s)
    s.set_defaults(fn=cmd_selftest)

    s = sub.add_parser("plot", help="figures to decide what to trust")
    s.add_argument("kind", choices=["selftest", "power", "coverage", "margins", "identity"],
                   help="selftest: can the rule see each gene (before a screen); power: are calls short of "
                        "reads (after a screen); coverage, margins, identity: why one run got its call")
    s.add_argument("--out", required=True, help="PNG or PDF; the values drawn go to a TSV of the same name")
    s.add_argument("--selftest", help="table from readloss selftest (plot selftest)")
    s.add_argument("--results", help="results table from test or batch (plot power, or with --run)")
    s.add_argument("--run", help="run accession or label in --results (per-run plots)")
    s.add_argument("--hits", help="a hit table, instead of --results and --run")
    s.add_argument("--label")
    s.add_argument("--genes"); s.add_argument("--proteins")
    s.add_argument("--family-fasta", help="the family FASTA, to map positions to one reference protein per gene")
    s.add_argument("--assignable")
    s.add_argument("--rule", choices=[assign.SYMMETRIC, assign.CENSUS], default=assign.SYMMETRIC)
    s.add_argument("--min-aln", type=int, default=None)
    s.add_argument("--margin", type=float, default=5.0)
    s.add_argument("--min-identity", type=float, default=None)
    s.add_argument("--nt", action="store_true")
    s.add_argument("--absent-max", type=float, default=0.10)
    s.add_argument("--present-min", type=float, default=0.30)
    s.add_argument("--min-expected", type=float, default=10.0)
    s.set_defaults(fn=cmd_plot)

    s = sub.add_parser("report", help="summary of a results table")
    s.add_argument("--results", required=True)
    s.add_argument("--out", help="summary TSV [stdout]")
    s.add_argument("--plot", help="PNG or PDF with the ratios (needs matplotlib)")
    s.add_argument("--absent-max", type=float, default=0.10)
    s.add_argument("--present-min", type=float, default=0.30)
    s.add_argument("--min-expected", type=float, default=10.0)
    s.set_defaults(fn=cmd_report)

    s = sub.add_parser("calibrate", help="sweep margin and identity on samples of known status")
    s.add_argument("--samples", required=True, help="TSV with columns label, hits, known_status (present/absent)")
    s.add_argument("--margins", default="0,2.5,5,7.5,10")
    s.add_argument("--identities", default="60,70,75,80,90")
    s.add_argument("--out")
    _common(s)
    s.set_defaults(fn=cmd_calibrate)

    s = sub.add_parser("sensitive", help="tBLASTn + BLASTx test (reproduces the genome tests)")
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument("--run", help="SRA run accession (streamed into the BLAST database)")
    g.add_argument("--reads", help="local FASTA of reads")
    g.add_argument("--blastx", help="existing BLASTx table of hit reads against the family")
    s.add_argument("--queries", help="candidate and control proteins for tBLASTn")
    s.add_argument("--family", help="FASTA with every protein of the family (BLASTx target)")
    s.add_argument("--spots", type=int, default=15_000_000)
    s.add_argument("--threads", type=int, default=4)
    s.add_argument("--workdir")
    s.add_argument("--keep-db", action="store_true", help="keep the BLAST database of the reads")
    s.add_argument("--label")
    s.add_argument("--out")
    _common(s)
    s.set_defaults(fn=cmd_sensitive)

    a = p.parse_args(argv)
    if a.cmd == "sensitive" and not a.blastx and not (a.queries and a.family):
        p.error("sensitive needs --queries and --family unless --blastx is given")
    a.fn(a)


if __name__ == "__main__":
    main()
