"""Detectability of each gene under the assignment rule (``readloss selftest``).

The ratio assumes that every read of the candidate and of the controls can be
assigned to its gene. That fails when family members are nearly identical, as
homeologs of an allopolyploid or recent duplicates: a read that does not cover a
position where they differ cannot be told apart and is discarded. If this happens
more to the candidate than to the controls, the ratio is biased, and if the
candidate can hardly ever be assigned, an absence call means nothing.

selftest measures this before any real run. It cuts every protein of the candidate
and of the controls into overlapping fragments the length of a read (50 aa for
150 nt reads), optionally adds substitutions to mimic a divergent species (one rate for all genes or
one per gene, since a gene on its way to loss often evolves faster than its controls, and
gamma-distributed rates across sites), searches
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


RATE_BLOCK = 20          # sites per block of shared rate (domains and linkers are contiguous)


def site_rates(n, alpha, rng, block=RATE_BLOCK):
    """Relative substitution rates of ``n`` sites, gamma distributed with mean 1 and shape ``alpha``.

    Real proteins mix conserved domains with fast, often disordered regions, and both are
    stretches of contiguous sites, so one rate is drawn per block of ``block`` sites; alpha
    around 0.5 to 1 is typical. ``alpha`` None or <= 0 gives equal rates.
    """
    if not alpha or alpha <= 0:
        return [1.0] * n
    out = []
    while len(out) < n:
        out += [rng.gammavariate(alpha, 1 / alpha)] * block
    return out[:n]


def mutate(seq, rate, rng, rates=None):
    if rate <= 0:
        return seq
    rates = rates or [1.0] * len(seq)
    return "".join(rng.choice(AA.replace(c, "")) if c in AA and rng.random() < min(1.0, rate * r) else c
                   for c, r in zip(seq, rates))


def fragments(family, seqs, length=50, step=5, divergence=0.0, seed=1, gene_divergence=None, alpha=None):
    """Fragments of every candidate and control protein, as (id, sequence); id is GENE::protein::start.

    The whole protein is mutated once, with ``gene_divergence[gene]`` (or ``divergence``) and
    gamma-distributed site rates, and then cut, so that fragments share the pattern of
    conserved and fast regions of the protein.
    """
    rng = random.Random(seed)
    tested = [family.candidate, *family.controls]
    gd = gene_divergence or {}
    out = []
    for pid, gene in family.protein_gene.items():
        if gene not in tested:
            continue
        if pid not in seqs:
            raise KeyError(f"protein {pid!r} of proteins.tsv is missing from the family FASTA")
        s = mutate(seqs[pid], gd.get(gene, divergence), rng, site_rates(len(seqs[pid]), alpha, rng))
        starts = range(0, max(1, len(s) - length + 1), step)
        for i in starts:
            out.append((f"{gene}{SEP}{pid}{SEP}{i + 1}", s[i:i + length]))
    return out


def read_divergence(path):
    """Table gene <tab> divergence (fraction of residues), with header."""
    import csv
    with open(path, newline="") as fh:
        return {r["gene"]: float(r["divergence"]) for r in csv.DictReader(fh, delimiter="\t")}


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
                threads=4, workdir=None, gene_divergence=None, alpha=None):
    """Fragment the proteins, search them with DIAMOND and evaluate the rule."""
    wd = workdir or tempfile.mkdtemp(prefix="readloss_selftest_")
    os.makedirs(wd, exist_ok=True)
    frags = fragments(family, read_fasta(family_faa), length, step, divergence, seed, gene_divergence, alpha)
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


NT = "ACGT"


def mutate_nt(seq, rate, rng):
    if rate <= 0:
        return seq
    return "".join(rng.choice(NT.replace(c, "")) if c in NT and rng.random() < rate else c for c in seq.upper())


def reads_nt(family, seqs, length=150, step=10, error=0.002, divergence=0.0, seed=1, gene_divergence=None):
    """Simulated reads of the candidate and control coding sequences, as (id, sequence).

    Each sequence is mutated once with its divergence (a divergent specimen), then cut into
    reads of ``length`` nt every ``step`` nt, on both strands, each with sequencing errors.
    """
    rng = random.Random(seed)
    tested = [family.candidate, *family.controls]
    gd = gene_divergence or {}
    comp = str.maketrans("ACGT", "TGCA")
    out = []
    for pid, gene in family.protein_gene.items():
        if gene not in tested:
            continue
        if pid not in seqs:
            raise KeyError(f"sequence {pid!r} of proteins.tsv is missing from the coding-sequence FASTA")
        s = mutate_nt(seqs[pid], gd.get(gene, divergence), rng)
        for i in range(0, max(1, len(s) - length + 1), step):
            r = mutate_nt(s[i:i + length], error, rng)
            if (i // step) % 2:
                r = r.translate(comp)[::-1]
            out.append((f"{gene}{SEP}{pid}{SEP}{i + 1}", r))
    return out


def run_nt(family, cds_fna, rule, length=150, step=10, error=0.002, divergence=0.0, seed=1, threads=4,
           workdir=None, gene_divergence=None):
    """Nucleotide selftest: simulated reads aligned with bwa-mem2 -a to the coding sequences."""
    wd = workdir or tempfile.mkdtemp(prefix="readloss_selftest_nt_")
    os.makedirs(wd, exist_ok=True)
    rs = reads_nt(family, read_fasta(cds_fna), length, step, error, divergence, seed, gene_divergence)
    q = os.path.join(wd, "reads.fa")
    with open(q, "w") as fh:
        for i, s in rs:
            fh.write(f">{i}\n{s}\n")
    db = pipeline.makedb_nt(cds_fna, os.path.join(wd, "family.fna"))
    tsv = pipeline.search_nt(q, db, os.path.join(wd, "reads.hits.tsv"), threads=threads)
    totals = Counter(origin(i) for i, _ in rs)
    return evaluate(assign.read_hits(tsv, family, rule), family, rule, totals)
