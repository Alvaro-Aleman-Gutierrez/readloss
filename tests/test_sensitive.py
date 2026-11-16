"""Sensitive mode reproduces Supplementary Table S2 from the original BLASTx tables."""
import csv
import gzip
import os
import random
import shutil

import pytest

from readloss import assign, family, sensitive, stats

D = os.path.join(os.path.dirname(__file__), "datos")
EX = os.path.join(os.path.dirname(__file__), "..", "ejemplos", "dcx")
SP = {"carbunculus": "Paedocypris carbunculus", "cerebrum_ctrl": "Danionella cerebrum",
      "dracula": "Danionella dracula", "furcatus": "Ictalurus furcatus",
      "misgurnus": "Misgurnus anguillicaudatus", "paedocypris": "Paedocypris micromegethes",
      "paramisgurnus": "Paramisgurnus dabryanus", "tinwini": "Danio tinwini"}
S2 = {r["species"]: r for r in csv.DictReader(open(os.path.join(D, "Table_S2.tsv")), delimiter="\t")}


@pytest.fixture(scope="module")
def fam():
    return family.load(os.path.join(EX, "genes_sensible.tsv"), os.path.join(EX, "proteins.tsv"))


@pytest.mark.parametrize("key", sorted(SP))
def test_table_s2(tmp_path, fam, key):
    p = tmp_path / "bx.tsv"
    with gzip.open(os.path.join(D, "sensible", f"{key}.blastx.tsv.gz"), "rt") as a, open(p, "w") as b:
        shutil.copyfileobj(a, b)
    s = assign.summarise(sensitive.assign_sensitive(p, fam), fam)
    r = S2[SP[key]]
    assert s["observed"] == int(r["reads_DCX"])
    for g in ("dclk1a", "dclk2a", "dclk3"):
        assert s[f"reads_{g}"] == int(r[f"reads_{g}"])
    assert round(s["expected"], 1) == pytest.approx(float(r["expected_DCX"]))
    assert round(s["ratio"], 3) == pytest.approx(float(r["obs_exp"]))
    assert s["status"] == r["DCX_status"]


def test_tie_goes_to_candidate(tmp_path, fam):
    p = tmp_path / "bx.tsv"
    p.write_text("r1\tDanrer_dclk2a_572548\t97\t39\t85.1\nr1\tHomsap_DCLK2_166614\t100\t39\t85.1\n")
    dcx = next(k for k, v in fam.protein_gene.items() if v == "DCX")
    p.write_text(p.read_text() + f"r1\t{dcx}\t90\t39\t85.1\n")
    c = sensitive.assign_sensitive(p, fam)
    assert c.by_gene["DCX"] == 1 and sum(c.by_gene.values()) == 1


def test_unique_ids():
    lines = [">a 1\n", "AC\n", ">a 1\n", "GT\n", ">b\n", "TT\n"]
    assert [l for l in sensitive.unique_ids(lines) if l.startswith(">")] == [">a\n", ">a_2\n", ">b\n"]


@pytest.mark.skipif(not all(shutil.which(t) for t in ("makeblastdb", "tblastn", "blastx", "blastdbcmd")),
                    reason="BLAST+ not installed")
def test_end_to_end_toy(tmp_path):
    rng = random.Random(7)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    cod = dict(zip(aa, ["GCT", "TGT", "GAT", "GAA", "TTT", "GGT", "CAT", "ATT", "AAA", "CTG",
                        "ATG", "AAT", "CCT", "CAA", "CGT", "TCT", "ACT", "GTT", "TGG", "TAT"]))
    L = {"cand": 350, "ctrlA": 700, "ctrlB": 800, "para": 400}
    P = {g: "".join(rng.choice(aa) for _ in range(n)) for g, n in L.items()}
    (tmp_path / "family.faa").write_text("".join(f">{g}_p1\n{s}\n" for g, s in P.items()))
    (tmp_path / "queries.faa").write_text("".join(f">{g}\n{P[g]}\n" for g in ("cand", "ctrlA", "ctrlB")))
    (tmp_path / "genes.tsv").write_text("gene\trole\tlength_aa\ncand\tcandidate\t350\n"
                                        "ctrlA\tcontrol\t700\nctrlB\tcontrol\t800\npara\tother\t0\n")
    (tmp_path / "proteins.tsv").write_text("protein_id\tgene\n" + "".join(f"{g}_p1\t{g}\n" for g in P))
    with open(tmp_path / "reads.fasta", "w") as fh:
        i = 0
        for g, s in P.items():
            for _ in range(len(s) // 5):
                st = rng.randrange(0, len(s) - 50)
                # same identifier twice, as fastq-dump writes the two mates of a spot
                seq = "".join(cod[c] for c in s[st:st + 50])
                fh.write(f">s{i} 1\n{seq}\n>s{i} 1\n{seq}\n")
                i += 1
    fam_ = family.load(tmp_path / "genes.tsv", tmp_path / "proteins.tsv")
    db = str(tmp_path / "reads_db")
    sensitive.build_reads_db(str(tmp_path / "reads.fasta"), db)
    bx = sensitive.run_search(db, str(tmp_path / "queries.faa"), str(tmp_path / "family.faa"), str(tmp_path))
    s = assign.summarise(sensitive.assign_sensitive(bx, fam_), fam_)
    assert s["observed"] == 2 * (350 // 5)            # both mates counted: the unit is the read
    assert s["status"] == stats.PRESENT and 0.9 < s["ratio"] < 1.1
