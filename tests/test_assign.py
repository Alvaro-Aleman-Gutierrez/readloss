"""The assignment gives the same counts as the original census script (parse_dmnd.py)."""
import os
import random
import re
import subprocess
import sys

import pytest

from readloss import assign, family

D = os.path.join(os.path.dirname(__file__), "datos")
# protein ids that the original script maps to each gene by name
PROT = {"DCX": ["Caraur_LOC113084161_a", "Danaes_DCX_b"],
        "dclk1a": ["Danrer_dclk1a_a", "Astmex_dclk1a_b"],
        "dclk2a": ["Danrer_dclk2a_a", "Lepocu_DCLK2A_b"],
        "other": ["Danrer_dclk3_a", "Danrer_dcdc2_a", "Danrer_dclk1b_a"]}


@pytest.fixture(scope="module")
def fam(tmp_path_factory):
    t = tmp_path_factory.mktemp("fam")
    (t / "genes.tsv").write_text("gene\trole\tlength_aa\nDCX\tcandidate\t357\n"
                                 "dclk1a\tcontrol\t770\ndclk2a\tcontrol\t826\nother\tother\t0\n")
    (t / "proteins.tsv").write_text("protein_id\tgene\n" +
                                    "".join(f"{p}\t{g}\n" for g, ps in PROT.items() for p in ps))
    return family.load(t / "genes.tsv", t / "proteins.tsv")


def synthetic(path, n_reads, seed):
    rng = random.Random(seed)
    allp = [p for ps in PROT.values() for p in ps]
    with open(path, "w") as fh:
        for i in range(n_reads):
            for _ in range(rng.randint(1, 5)):
                s = rng.choice(allp)
                bs = rng.choice([rng.uniform(20, 200), 100.0, 103.0, 105.0])   # force ties and near-margins
                fh.write(f"read{i}\t{s}\t{rng.uniform(40, 100):.1f}\t{rng.randint(20, 60)}\t{bs:.1f}\n")


def original(path):
    out = subprocess.run([sys.executable, os.path.join(D, "parse_dmnd_original.py"), str(path), "x"],
                         capture_output=True, text=True, check=True).stdout
    return {k: v for k, v in re.findall(r"(\w+)=([^\t\n]+)", out)}


@pytest.mark.parametrize("seed", range(20))
def test_same_as_original(tmp_path, fam, seed):
    f = tmp_path / "hits.tsv"
    synthetic(f, 400, seed)
    o = original(f)
    rule = assign.Rule(mode=assign.CENSUS)
    c = assign.assign(assign.read_hits(f, fam, rule), fam, rule)
    s = assign.summarise(c, fam)
    assert s["observed"] == int(o["DCX"])
    assert s["reads_dclk1a"] == int(o["dclk1a"])
    assert s["reads_dclk2a"] == int(o["dclk2a"])
    assert s["candidate_discarded"] == int(o["DCX_descartadas"])
    assert f"{s['expected']:.1f}" == o["esperado"]


def test_subsample_is_thinning(tmp_path, fam):
    f = tmp_path / "hits.tsv"
    synthetic(f, 4000, 1)
    hits = assign.read_hits(f, fam)
    keep = assign.subsample_ids(hits, 0.25, seed=3)
    assert 0.22 < len(keep) / len(hits) < 0.28
    assert keep == assign.subsample_ids(hits, 0.25, seed=3)


def test_unknown_protein(tmp_path, fam):
    f = tmp_path / "hits.tsv"
    f.write_text("r1\tUnknown_prot\t90\t40\t80\n")
    with pytest.raises(KeyError):
        assign.read_hits(f, fam)


def test_subsample_replicates_independent(tmp_path, fam):
    f = tmp_path / "hits.tsv"
    synthetic(f, 4000, 2)
    hits = assign.read_hits(f, fam)
    a = assign.subsample_ids(hits, 0.5, seed=1)
    b = assign.subsample_ids(hits, 0.5, seed=2)
    overlap = len(a & b) / len(hits)
    assert 0.2 < overlap < 0.3          # independent halves share about a quarter of reads
