"""Sensitive mode: tBLASTn of candidate and control proteins against the reads, then BLASTx
of every read with a hit against the whole family.

This is the version used for the genomes of the doublecortin study (Supplementary Table S2).
It is slower than the DIAMOND mode but finds more divergent reads, so it suits distant
queries or runs with few reads.

Differences from the DIAMOND mode, kept so that the published counts are reproduced:
  * the unit is the read: the two mates of a spot get distinct identifiers;
  * every alignment counts, with no minimum length, bit margin or identity filter;
  * a read goes to its best BLASTx hit, and a tie between the candidate and another member
    is counted for the candidate, which is conservative for a test of absence.
"""
from collections import Counter
import os
import re
import subprocess

from . import assign
from .pipeline import which

TBLASTN_FIELDS = "6 qseqid sseqid pident length"
BLASTX_FIELDS = "6 qseqid sseqid pident length bitscore sstart send"


def unique_ids(lines):
    """Yield FASTA lines with repeated identifiers made unique (id, id_2, id_3 ...)."""
    seen = Counter()
    for line in lines:
        if line.startswith(">"):
            head = line[1:].rstrip("\n")
            ident, _, rest = head.partition(" ")
            seen[ident] += 1
            if seen[ident] > 1:
                ident = f"{ident}_{seen[ident]}"
            line = f">{ident}\n"
        yield line


def build_reads_db(source, out_db, run=None, spots=None, log=None):
    """BLAST nucleotide database of the reads, from a FASTA file or streamed from the SRA.

    Identifiers are made unique on the way, so that hit reads can be pulled back with
    blastdbcmd.
    """
    mk = which("makeblastdb")
    cmd = [mk, "-in", "-", "-dbtype", "nucl", "-parse_seqids", "-out", out_db, "-title", out_db]
    mkp = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, text=True)
    written = None
    if run:
        fq = which("fastq-dump", "READLOSS_FASTQ_DUMP")
        dump = subprocess.Popen([fq, "--fasta", "0", "--split-spot", "--skip-technical", "-X", str(spots),
                                 "--stdout", run], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for line in unique_ids(dump.stdout):
            mkp.stdin.write(line)
        err = dump.stderr.read()
        dump.wait()
        if log:
            with open(log, "a") as fh:
                fh.write(err)
        m = re.search(r"^Written (\d+) spots", err, re.M)
        written = int(m.group(1)) if m else 0
        if dump.returncode != 0 or not written:
            mkp.kill()
            raise RuntimeError(f"download of {run} failed")
    else:
        with open(source) as fh:
            for line in unique_ids(fh):
                mkp.stdin.write(line)
    mkp.stdin.close()
    if mkp.wait() != 0:
        raise RuntimeError("makeblastdb failed")
    return written


def run_search(reads_db, queries_faa, family_faa, workdir, threads=4):
    """tBLASTn of the queries against the reads, then BLASTx of the hit reads against the family."""
    tbn = os.path.join(workdir, "tblastn.tsv")
    ids = os.path.join(workdir, "hit_ids.txt")
    reads = os.path.join(workdir, "hit_reads.fna")
    famdb = os.path.join(workdir, "family_prot")
    bx = os.path.join(workdir, "blastx.tsv")
    subprocess.run([which("tblastn"), "-query", queries_faa, "-db", reads_db, "-evalue", "1e-5",
                    "-outfmt", TBLASTN_FIELDS, "-max_target_seqs", "100000", "-seg", "no",
                    "-num_threads", str(threads), "-out", tbn], check=True)
    hit = sorted({l.split("\t")[1] for l in open(tbn) if l.strip()})
    with open(ids, "w") as fh:
        fh.write("".join(h + "\n" for h in hit))
    if not hit:
        open(bx, "w").close()
        return bx
    subprocess.run([which("blastdbcmd"), "-db", reads_db, "-entry_batch", ids, "-outfmt", "%f",
                    "-out", reads], check=True)
    subprocess.run([which("makeblastdb"), "-in", family_faa, "-dbtype", "prot", "-out", famdb],
                   check=True, stdout=subprocess.DEVNULL)
    subprocess.run([which("blastx"), "-query", reads, "-db", famdb, "-evalue", "1e-5",
                    "-outfmt", BLASTX_FIELDS, "-max_target_seqs", "5", "-max_hsps", "1",
                    "-num_threads", str(threads), "-out", bx], check=True)
    return bx


def assign_sensitive(blastx_tsv, family):
    """Best BLASTx hit per read; ties that include the candidate go to the candidate."""
    best = {}
    cand = family.candidate
    for line in open(blastx_tsv):
        if not line.strip():
            continue
        f = line.split()
        q, s, pid, ln, bs = f[:5]
        span = (min(int(f[5]), int(f[6])), max(int(f[5]), int(f[6]))) if len(f) >= 7 else None
        bs, g = float(bs), family.gene_of(s)
        if q not in best or bs > best[q][0]:
            best[q] = (bs, g, float(pid), span)
        elif bs == best[q][0] and g == cand:
            best[q] = (bs, g, float(pid), span)
    c = assign.Counts()
    for bs, g, pid, span in best.values():
        c.by_gene[g] += 1
        c.identity[g].append(pid)
        if g == cand:
            c.candidate_identity.append(pid)
            if span:
                c.candidate_spans.append(span)
    return c
