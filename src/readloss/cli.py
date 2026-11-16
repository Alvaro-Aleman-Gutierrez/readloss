"""Command line interface.

readloss makedb     build the DIAMOND database of the gene family
readloss test       test one run (SRA accession, local FASTA or existing hit table)
readloss batch      test many runs, resumable, with automatic deeper repeat of intermediates
readloss subsample  ratio and limits against read depth, from one hit table
readloss calibrate  sweep the assignment rule on samples of known status
readloss sensitive  slower tBLASTn + BLASTx test for divergent queries or few reads
readloss depth      reads needed, from a short pilot run or an existing hit table
"""
import argparse
import csv
import os
import sys
import tempfile

from . import __version__, assign, family as fam, pipeline, sensitive, stats

FMT = {"expected": "{:.1f}", "ratio": "{:.3f}", "ci_low": "{:.3f}", "ci_high": "{:.3f}",
       "p_le_obs_if_present": "{:.2g}", "mean_identity_candidate": "{:.1f}"}


def _fmt(row):
    return {k: (FMT[k].format(v) if k in FMT and isinstance(v, float) else v) for k, v in row.items()}


def _write(rows, path):
    """Write rows as TSV to stdout or append to ``path`` (header only for a new file)."""
    rows = [_fmt(r) for r in rows]
    if not rows:
        return
    new = not path or not os.path.exists(path) or os.path.getsize(path) == 0
    out = open(path, "a", newline="") if path else sys.stdout
    w = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore")
    if new:
        w.writeheader()
    w.writerows(rows)
    if path:
        out.close()


def _common(p):
    p.add_argument("--genes", required=True, help="genes.tsv (gene, role, length_aa)")
    p.add_argument("--proteins", required=True, help="proteins.tsv (protein_id, gene)")
    p.add_argument("--min-aln", type=int, default=30, help="minimum alignment length in aa [30]")
    p.add_argument("--margin", type=float, default=5.0, help="bits over the best non-candidate hit [5]")
    p.add_argument("--min-identity", type=float, default=75.0, help="identity for candidate reads, percent [75]")
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


def _setup(a):
    f = fam.load(a.genes, a.proteins)
    rule = assign.Rule(a.min_aln, a.margin, a.min_identity)
    th = stats.Thresholds(a.absent_max, a.present_min, a.min_expected)
    return f, rule, th


def _from_hits(hits_tsv, f, rule, th):
    hits = assign.read_hits(hits_tsv, f, rule)
    return assign.summarise(assign.assign(hits, f, rule), f, th)


def _one_pass(run, spots, a, f, rule, th, wd):
    fasta = os.path.join(wd, f"{run}.{spots}.fasta")
    tsv = os.path.join(wd, f"{run}.{spots}.hits.tsv")
    spots_file = tsv + ".spots"
    if not os.path.exists(tsv):
        log = os.path.join(wd, "download.log")
        if a.stream:
            got = pipeline.fetch_search_stream(run, spots, a.db, tsv, threads=a.threads, log=log,
                                               allow_short=a.allow_short)
        else:
            got = pipeline.fetch(run, spots, fasta, log=log, allow_short=a.allow_short)
            pipeline.search(fasta, a.db, tsv, threads=a.threads)
            if not a.keep_reads:
                os.remove(fasta)
        with open(spots_file, "w") as fh:
            fh.write(f"{got}\n")
    else:
        got = int(open(spots_file).read()) if os.path.exists(spots_file) else spots
    row = _from_hits(tsv, f, rule, th)
    row["spots"] = got
    row["hits_table"] = tsv
    return row


def _test_run(run, label, a, f, rule, th, wd):
    first = _one_pass(run, a.spots, a, f, rule, th, wd)
    final, repeat = first, False
    if first["status"] == stats.INTERMEDIATE and a.spots_max and a.spots_max > a.spots:
        final, repeat = _one_pass(run, a.spots_max, a, f, rule, th, wd), True
    row = {"label": label, "run": run}
    row.update(final)
    row["repeated_deeper"] = repeat
    if repeat:
        row.update({"first_spots": first["spots"], "first_observed": first["observed"],
                    "first_expected": first["expected"], "first_ratio": first["ratio"]})
    return row


def cmd_makedb(a):
    print(pipeline.makedb(a.family, a.out))


def cmd_test(a):
    f, rule, th = _setup(a)
    if a.hits:
        row = {"label": a.label or os.path.basename(a.hits)}
        row.update(_from_hits(a.hits, f, rule, th))
    elif a.reads:
        if not a.db:
            sys.exit("--db is required with --reads")
        wd = a.workdir or tempfile.mkdtemp(prefix="readloss_")
        tsv = os.path.join(wd, os.path.basename(a.reads) + ".hits.tsv")
        pipeline.search(a.reads, a.db, tsv, threads=a.threads)
        row = {"label": a.label or os.path.basename(a.reads)}
        row.update(_from_hits(tsv, f, rule, th))
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
            row = _test_run(r["run"], r.get("label") or r["run"], a, f, rule, th, wd)
        except Exception as e:                       # keep going, log the failure
            with open(os.path.join(wd, "failed.tsv"), "a") as fh:
                fh.write(f"{r['run']}\t{e}\n")
            continue
        _write([row], a.out)


def cmd_subsample(a):
    f, rule, th = _setup(a)
    hits = assign.read_hits(a.hits, f, rule)
    rows = []
    for fr in [float(x) for x in a.fractions.split(",")]:
        for rep in range(1, a.replicates + 1):
            keep = assign.subsample_ids(hits, fr, seed=rep)
            row = {"fraction": fr, "replicate": rep}
            if a.spots:
                row["spots"] = round(a.spots * fr)
            row.update(assign.summarise(assign.assign(hits, f, rule, reads=keep), f, th))
            rows.append(row)
    _write(rows, a.out)


def cmd_depth(a):
    """Expected candidate reads per million spots and the depth that reaches the target."""
    f, rule, th = _setup(a)
    if a.hits:
        if not a.spots:
            sys.exit("--spots (spots behind the hit table) is required with --hits")
        row = _from_hits(a.hits, f, rule, th)
        spots, label = a.spots, a.label or os.path.basename(a.hits)
    else:
        if not a.db:
            sys.exit("--db is required with --run")
        wd = a.workdir or tempfile.mkdtemp(prefix="readloss_")
        os.makedirs(wd, exist_ok=True)
        a.stream = True
        row = _one_pass(a.run, a.pilot_spots, a, f, rule, th, wd)
        spots, label = row["spots"], a.label or a.run
    rate = row["expected"] / spots * 1e6
    need = a.safety * th.min_expected / rate * 1e6 if rate > 0 else float("nan")
    out = {"label": label, "spots_used": spots, "expected": row["expected"],
           "expected_per_million_spots": round(rate, 3), "target_expected": th.min_expected,
           "safety": a.safety, "recommended_spots": int(-(-need // 1e6) * 1e6) if need == need else "nan"}
    _write([out], a.out)


def cmd_calibrate(a):
    f, _, th = _setup(a)
    with open(a.samples) as fh:
        samples = [r for r in csv.DictReader(fh, delimiter="\t")]
    rows = []
    for m in [float(x) for x in a.margins.split(",")]:
        for i in [float(x) for x in a.identities.split(",")]:
            rule = assign.Rule(a.min_aln, m, i)
            ratios = {"present": [], "absent": []}
            for s in samples:
                res = _from_hits(s["hits"], f, rule, th)
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


def cmd_sensitive(a):
    f = fam.load(a.genes, a.proteins)
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
    row.update(assign.summarise(sensitive.assign_sensitive(bx, f), f, th))
    _write([row], a.out)


def main(argv=None):
    p = argparse.ArgumentParser(prog="readloss", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version=f"readloss {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("makedb", help="build the DIAMOND database of the family")
    s.add_argument("--family", required=True, help="FASTA with every protein of the family")
    s.add_argument("--out", required=True, help="output database name")
    s.set_defaults(fn=cmd_makedb)

    s = sub.add_parser("test", help="test one run")
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument("--run", help="SRA run accession")
    g.add_argument("--reads", help="local FASTA of reads")
    g.add_argument("--hits", help="existing DIAMOND table (qseqid sseqid pident length bitscore)")
    s.add_argument("--label")
    s.add_argument("--out", help="output TSV [stdout]")
    _common(s); _search_opts(s)
    s.set_defaults(fn=cmd_test)

    s = sub.add_parser("batch", help="test many runs")
    s.add_argument("--runs", required=True, help="TSV with columns run and label")
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
    s.add_argument("--allow-short", action="store_true"); s.add_argument("--label"); s.add_argument("--out")
    _common(s)
    s.set_defaults(fn=cmd_depth, keep_reads=False)

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
