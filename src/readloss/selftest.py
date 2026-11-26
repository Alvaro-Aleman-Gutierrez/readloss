"""Detectability of each gene under the assignment rule (``readloss selftest``).

The ratio assumes that every read of the candidate and of the controls can be
assigned to its gene. That fails when family members are nearly identical, as
homeologs of an allopolyploid or recent duplicates: a read that does not cover a
position where they differ cannot be told apart and is discarded. If this happens
more to the candidate than to the controls, the ratio is biased, and if the
candidate can hardly ever be assigned, an absence call means nothing.

selftest measures this before any real run. It cuts every protein of the candidate
and of the controls into overlapping fragments the length of a read (50 aa for
150 nt reads), optionally adds substitutions to mimic a divergent species, searches
them against the whole family with the settings of the real search and applies the
same rule. For each gene it reports the fraction of fragments assigned back to it.
That fraction is the assignable fraction used by ``--assignable`` to correct the
expected count:

    expected = sum(control reads) / sum(control length x assignable) x candidate length x assignable

Reads simulated in another way, for example at the nucleotide level, can be given
instead as a hit table whose read identifiers start with the gene of origin,
``GENE::anything``.
"""
import os
import random
import tempfile
from collections import Counter

from . import assign, pipeline

AA = "ACDEFGHIKLMNPQRSTVWY"
SEP = "::"


def read_fasta(path):
    seqs, k = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                k = line[1:].split()[0]
                seqs[k] = []
            elif k is not None:
                seqs[k].append(line)
    return {a: "".join(b) for a, b in seqs.items()}


def mutate(seq, rate, rng):
    if rate <= 0:
        return seq
    return "".join(rng.choice(AA.replace(c, "")) if c in AA and rng.random() < rate else c for c in seq)


def fragments(family, seqs, length=50, step=5, divergence=0.0, seed=1):
    """Fragments of every candidate and control protein, as (id, sequence); id is GENE::protein::start."""
    rng = random.Random(seed)
    tested = [family.candidate, *family.controls]
    out = []
    for pid, gene in family.protein_gene.items():
        if gene not in tested:
            continue
        if pid not in seqs:
            raise KeyError(f"protein {pid!r} of proteins.tsv is missing from the family FASTA")
        s = seqs[pid]
        starts = range(0, max(1, len(s) - length + 1), step)
        for i in starts:
            frag = mutate(s[i:i + length], divergence, rng)
            out.append((f"{gene}{SEP}{pid}{SEP}{i + 1}", frag))
    return out


def origin(read_id):
    return read_id.split(SEP, 1)[0]


def evaluate(hits, family, rule, totals):
    """Per tested gene: where its fragments went under the rule.

    ``hits`` is the output of assign.read_hits; ``totals`` maps gene -> fragments made,
    including those without any hit, which count as not assignable.
    """
    by_origin = {}
    for q in hits:
        by_origin.setdefault(origin(q), set()).add(q)
    rows = []
    for g in [family.candidate, *family.controls]:
        ids = by_origin.get(g, set())
        c = assign.assign(hits, family, rule, reads=ids) if ids else assign.Counts()
        n = totals.get(g, len(ids))
        own = c.by_gene.get(g, 0)
        to_cand = c.by_gene.get(family.candidate, 0) if g != family.candidate else own
        to_ctrl = sum(c.by_gene.get(x, 0) for x in family.controls if x != g)
        counted = sum(v for k, v in c.by_gene.items() if k != assign.DISCARDED)
        discarded = len(ids) - counted
        rows.append({
            "gene": g,
            "role": "candidate" if g == family.candidate else "control",
            "fragments": n,
            "with_hits": len(ids),
            "assigned_self": own,
            "assigned_candidate": to_cand if g != family.candidate else 0,
            "assigned_other_controls": to_ctrl,
            "assigned_other_genes": counted - own - (to_cand if g != family.candidate else 0) - to_ctrl,
            "discarded": discarded,
            "no_hit": n - len(ids),
            "assignable_fraction": own / n if n else 0.0,
        })
    return rows


def run_protein(family, family_faa, rule, db=None, length=50, step=5, divergence=0.0, seed=1,
                threads=4, workdir=None):
    """Fragment the proteins, search them with DIAMOND and evaluate the rule."""
    wd = workdir or tempfile.mkdtemp(prefix="readloss_selftest_")
    os.makedirs(wd, exist_ok=True)
    frags = fragments(family, read_fasta(family_faa), length, step, divergence, seed)
    q = os.path.join(wd, "fragments.faa")
    with open(q, "w") as fh:
        for i, s in frags:
            fh.write(f">{i}\n{s}\n")
    if not db:
        db = pipeline.makedb(family_faa, os.path.join(wd, "family"))
    tsv = pipeline.search_proteins(q, db, os.path.join(wd, "fragments.hits.tsv"), threads=threads)
    totals = Counter(origin(i) for i, _ in frags)
    return evaluate(assign.read_hits(tsv, family, rule), family, rule, totals)


def run_hits(family, hits_tsv, rule, totals=None):
    """Evaluate the rule on an existing hit table of simulated reads named GENE::..."""
    return evaluate(assign.read_hits(hits_tsv, family, rule), family, rule, totals or {})
