"""The Snakemake workflow runs end to end on the example hit tables (skipped without Snakemake)."""
import csv
import gzip
import os
import re
import shutil
import subprocess

import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")
EX = os.path.join(ROOT, "ejemplos", "dcx")


@pytest.mark.skipif(not shutil.which("snakemake"), reason="Snakemake not installed")
def test_workflow(tmp_path):
    hits = tmp_path / "hits"
    hits.mkdir()
    for f in os.listdir(os.path.join(EX, "hits")):
        run = re.search(r"([A-Z]RR\d+)\.hits\.tsv\.gz$", f).group(1)
        with gzip.open(os.path.join(EX, "hits", f), "rt") as a, open(hits / f"{run}.hits.tsv", "w") as b:
            shutil.copyfileobj(a, b)
    cfg = tmp_path / "config.yaml"
    cfg.write_text(f"runs: {EX}/runs.tsv\nfamily: {EX}/family.faa\ngenes: {EX}/genes.tsv\n"
                   f"proteins: {EX}/proteins.tsv\noutdir: {tmp_path}/out\nhits_dir: {hits}\nthreads: 1\n")
    subprocess.run(["snakemake", "-s", os.path.join(ROOT, "workflow", "Snakefile"), "--configfile", str(cfg),
                    "--cores", "2", "-q", "--directory", str(tmp_path)], check=True)
    rows = list(csv.DictReader(open(tmp_path / "out" / "results.tsv"), delimiter="\t"))
    status = {r["label"]: r["status"] for r in rows}
    assert status["Paedocypris_micromegethes_absent"] == "absent"
    assert status["Danionella_cerebrum_present"] == "present"
    assert (tmp_path / "out" / "ratios.pdf").stat().st_size > 0
