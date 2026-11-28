"""External steps: partial download from the SRA and DIAMOND search."""
import os
import re
import shutil
import subprocess
import time

DIAMOND_FIELDS = ["qseqid", "sseqid", "pident", "length", "bitscore", "sstart", "send",
                  "qseq_translated", "full_qseq"]
# the last two give the translated aligned read (stop codons) and the whole read (GC, length)

RUNINFO_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&rettype=runinfo&tool=readloss&id="
WGS_STRATEGIES = {"WGS", "WGA"}
RANDOM_SELECTIONS = {"RANDOM", "RANDOM PCR", "SIZE FRACTIONATION", "UNSPECIFIED", "OTHER"}
GENOMIC_SOURCES = {"GENOMIC", "GENOMIC SINGLE CELL"}


def parse_runinfo(text, run):
    """The row of ``run`` in an SRA runinfo answer, XML or CSV, as a dict (empty if absent)."""
    t = text.lstrip()
    if t.startswith("<"):
        import xml.etree.ElementTree as ET
        root = ET.fromstring(t)
        for row in root.iter("Row"):
            d = {c.tag: (c.text or "") for c in row}
            if d.get("Run") == run:
                return d
        return {}
    import csv
    import io
    for r in csv.DictReader(io.StringIO(text)):
        if r.get("Run") == run:
            return r
    return {}


def run_info(run, timeout=30, retries=3, wait=5):
    """SRA runinfo of a run (library strategy, selection, source, layout, spots, read length, sex)."""
    import urllib.request
    last = None
    for i in range(retries):
        try:
            with urllib.request.urlopen(RUNINFO_URL + run, timeout=timeout) as r:
                return parse_runinfo(r.read().decode("utf-8", "replace"), run)
        except Exception as e:                   # network hiccup: retry, then give up quietly
            last = e
            time.sleep(wait)
    raise RuntimeError(f"SRA metadata of {run} could not be read: {last}")


EXPERIMENT_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&retmode=xml&tool=readloss&id="
_PCR_FREE = re.compile(r"pcr[\s_-]*free|amplification[\s_-]*free|\bno[\s_-]*pcr\b|without\s+pcr", re.I)
_PCR = re.compile(r"\bpcr\b|nextera|truseq\s+nano|nebnext\s+ultra|kapa\s+hyper|amplified|amplification", re.I)


def pcr_free(text):
    """yes, no or unknown, from the library description of a run.

    Libraries made without PCR cover GC-rich DNA far more evenly, so the GC diagnostics
    matter less for them.
    """
    if not text:
        return "unknown"
    if _PCR_FREE.search(text):
        return "yes"
    if _PCR.search(text):
        return "no"
    return "unknown"


def parse_experiment(text):
    """Library description fields of an SRA experiment XML, joined (empty if none)."""
    import xml.etree.ElementTree as ET
    try:
        root = ET.fromstring(text.lstrip())
    except ET.ParseError:
        return ""
    parts = []
    for tag in ("DESIGN_DESCRIPTION", "LIBRARY_CONSTRUCTION_PROTOCOL", "LIBRARY_NAME"):
        for e in root.iter(tag):
            if e.text and e.text.strip():
                parts.append(e.text.strip())
    return " | ".join(parts)


def parse_aligned(text, run=None):
    """True if an SRA package XML describes ``run`` as submitted aligned (BAM or CRAM, with AlignInfo).

    One record can hold every run of an experiment, so only the RUN element of ``run`` is read
    (all of them when ``run`` is None).
    """
    import xml.etree.ElementTree as ET
    try:
        root = ET.fromstring(text.lstrip())
    except ET.ParseError:
        return False
    runs = [e for e in root.iter("RUN") if run is None or e.get("accession") == run]
    for e in runs or ([root] if run is None else []):
        if any(True for _ in e.iter("AlignInfo")):
            return True
        if any((f.get("semantic_name") or "").lower() in ("bam", "cram") and f.get("supertype") == "Original"
               for f in e.iter("SRAFile")):
            return True
    return False


def library_package(run, timeout=30):
    """(library description, submitted aligned) of a run, from one SRA package record."""
    import urllib.request
    with urllib.request.urlopen(EXPERIMENT_URL + run, timeout=timeout) as r:
        text = r.read().decode("utf-8", "replace")
    return parse_experiment(text), parse_aligned(text, run)


def library_description(run, timeout=30):
    """Design and library-construction text of a run, from the SRA experiment record."""
    import urllib.request
    with urllib.request.urlopen(EXPERIMENT_URL + run, timeout=timeout) as r:
        return parse_experiment(r.read().decode("utf-8", "replace"))


def library_issues(info):
    """Reasons why a run is not a random whole-genome library (empty list if it is)."""
    out = []
    st = (info.get("LibraryStrategy") or "").upper()
    se = (info.get("LibrarySelection") or "").upper()
    so = (info.get("LibrarySource") or "").upper()
    if st and st not in WGS_STRATEGIES:
        out.append(f"strategy {st}")
    if se and se not in RANDOM_SELECTIONS:
        out.append(f"selection {se}")
    if so and so not in GENOMIC_SOURCES:
        out.append(f"source {so}")
    return out


def mate_length(info):
    """Mean length of one read (nt) from runinfo; avgLength sums both mates of a pair."""
    try:
        L = float(info.get("avgLength") or "nan")
    except ValueError:
        return float("nan")
    return L / 2 if (info.get("LibraryLayout") or "").upper() == "PAIRED" else L


def spot_ranges(spots, total, blocks):
    """Spot ranges spread along a run: ``blocks`` blocks that together hold ``spots`` spots.

    The first spots of a run can all come from one lane or one library of several merged
    in a run, and the edges of a flow cell give poorer reads, so spreading the sample
    along the run makes it closer to a random sample.
    """
    if blocks <= 1 or not total or total <= spots:
        return [(1, spots)]
    per = spots // blocks
    step = total // blocks
    rs = [(1 + i * step, i * step + per) for i in range(blocks)]
    last = rs[-1]
    rs[-1] = (last[0], last[1] + spots - per * blocks)     # remainder in the last block
    return rs


def _dump_args(run, ranges):
    if ranges is None:
        return None
    return [["-N", str(a), "-X", str(b)] for a, b in ranges]


def _written(err):
    return sum(int(x) for x in re.findall(r"^Written (\d+) spots", err, re.M))


def which(tool, env=None):
    path = os.environ.get(env) if env else None
    path = path or shutil.which(tool)
    if not path:
        raise FileNotFoundError(f"{tool} not found in PATH" + (f" (or set {env})" if env else ""))
    return path


def fetch(run, spots, out_fasta, retries=3, wait=30, log=None, allow_short=False, ranges=None):
    """Download the first ``spots`` spots of a run as FASTA, mates split.

    The download is accepted only if fastq-dump reports exactly ``spots`` written,
    so a truncated transfer is never analysed. With ``allow_short`` a run that is
    smaller than requested is accepted as it is (check the log, because a broken
    transfer that exits cleanly would look the same).
    """
    exe = which("fastq-dump", "READLOSS_FASTQ_DUMP")
    base = [exe, "--fasta", "0", "--split-spot", "--skip-technical"]
    cmds = [base + ["-X", str(spots), "--stdout", run]] if not ranges else \
        [base + r + ["--stdout", run] for r in _dump_args(run, ranges)]
    for i in range(1, retries + 1):
        err, rc = "", 0
        with open(out_fasta, "w") as out:
            for cmd in cmds:
                p = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
                err += p.stderr
                rc = rc or p.returncode
        if log:
            with open(log, "a") as fh:
                fh.write(err)
        n = _written(err)
        if rc == 0 and n == spots:
            return n
        if allow_short and rc == 0 and 0 < n < spots:
            return n
        os.remove(out_fasta)
        if i < retries:
            time.sleep(wait)
    raise RuntimeError(f"download of {run} failed after {retries} attempts")


def fetch_search_stream(run, spots, db, out_tsv, threads=4, retries=3, wait=30, log=None,
                        allow_short=False, ranges=None, search_fn=None):
    """Download and search at once through a named pipe, so no reads are written to disk.

    The hit table is kept only if fastq-dump reports the full number of spots.
    """
    fq = which("fastq-dump", "READLOSS_FASTQ_DUMP")
    work = os.path.dirname(os.path.abspath(out_tsv))
    for i in range(1, retries + 1):
        fifo = os.path.join(work, f".{run}.{os.getpid()}.fifo")
        if os.path.exists(fifo):
            os.remove(fifo)
        os.mkfifo(fifo)
        if ranges:
            body = "; ".join(f'"$0" --fasta 0 --split-spot --skip-technical -N {lo} -X {hi} --stdout {run} '
                             f'|| exit 1' for lo, hi in ranges)
            script = f'{{ {body}; }} > "$1"'
        else:
            script = f'exec "$0" --fasta 0 --split-spot --skip-technical -X {spots} --stdout {run} > "$1"'
        dump = subprocess.Popen(["sh", "-c", script, fq, fifo], stderr=subprocess.PIPE, text=True)
        try:
            (search_fn or search)(fifo, db, out_tsv, threads=threads)
            ok_search = True
        except subprocess.CalledProcessError:
            ok_search = False
            dump.kill()                      # nobody reads the pipe any more
        finally:
            err = dump.communicate()[1]
            os.remove(fifo)
        if log:
            with open(log, "a") as fh:
                fh.write(err)
        n = _written(err)
        if ok_search and dump.returncode == 0 and (n == spots or (allow_short and 0 < n < spots)):
            return n
        if os.path.exists(out_tsv):
            os.remove(out_tsv)
        if i < retries:
            time.sleep(wait)
    raise RuntimeError(f"streamed download and search of {run} failed after {retries} attempts")


def makedb(family_faa, out_db, panel_faa=None):
    """DIAMOND database of the family, optionally with a panel of single-copy genes.

    With a panel, the combined FASTA and ``<out>.panel.tsv`` (protein_id, gene, length_aa)
    are written next to the database; pass the latter to ``--panel``.
    """
    exe = which("diamond", "READLOSS_DIAMOND")
    src = family_faa
    if panel_faa:
        stem = out_db[:-5] if out_db.endswith(".dmnd") else out_db
        src = stem + ".family_panel.faa"
        fam_ids = set()
        with open(family_faa) as fh:
            fam_ids = {l[1:].split()[0] for l in fh if l.startswith(">")}
        rows, seq, pid = [], [], None
        with open(src, "w") as out, open(family_faa) as fh:
            out.write(fh.read().rstrip("\n") + "\n")
            with open(panel_faa) as ph:
                for line in list(ph) + [">"]:
                    if line.startswith(">"):
                        if pid is not None:
                            s_ = "".join(seq)
                            if pid in fam_ids:
                                raise ValueError(f"panel protein {pid!r} is also in the family FASTA")
                            out.write(f">{pid}\n{s_}\n")
                            rows.append((pid, f"panel:{pid}", len(s_)))
                        pid, seq = (line[1:].split()[0] if len(line) > 2 else None), []
                    else:
                        seq.append(line.strip())
        with open(stem + ".panel.tsv", "w") as t:
            t.write("protein_id\tgene\tlength_aa\n" + "".join(f"{a}\t{g}\t{n}\n" for a, g, n in rows))
    subprocess.run([exe, "makedb", "--in", src, "-d", out_db, "--quiet"], check=True)
    return out_db if out_db.endswith(".dmnd") else out_db + ".dmnd"


_SEQ_FIELDS = {}

_PROBE_PROT = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ"
_PROBE_CODON = {"M": "ATG", "K": "AAA", "T": "ACC", "A": "GCT", "Y": "TAT", "I": "ATT", "Q": "CAA", "R": "CGT",
                "S": "TCT", "F": "TTT", "V": "GTT", "H": "CAT", "L": "CTG", "E": "GAA", "G": "GGT", "P": "CCT",
                "D": "GAT", "N": "AAT", "W": "TGG", "C": "TGT"}


def seq_fields_supported(exe):
    """Whether this DIAMOND accepts qseq_translated and full_qseq (probed once with a tiny search)."""
    if exe in _SEQ_FIELDS:
        return _SEQ_FIELDS[exe]
    import tempfile
    ok = False
    with tempfile.TemporaryDirectory(prefix="readloss_probe_") as d:
        faa, fa, db = os.path.join(d, "p.faa"), os.path.join(d, "r.fa"), os.path.join(d, "p")
        with open(faa, "w") as fh:
            fh.write(f">P\n{_PROBE_PROT}\n")
        with open(fa, "w") as fh:
            fh.write(">r\n" + "".join(_PROBE_CODON[c] for c in _PROBE_PROT[:50]) + "\n")
        try:
            subprocess.run([exe, "makedb", "--in", faa, "-d", db, "--quiet"], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            r = subprocess.run([exe, "blastx", "-q", fa, "-d", db, "--quiet", "-f", "6", *DIAMOND_FIELDS,
                                "-o", os.path.join(d, "o.tsv")], stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL)
            ok = r.returncode == 0
        except (OSError, subprocess.CalledProcessError):
            ok = False
    _SEQ_FIELDS[exe] = ok
    return ok


def search(reads_fasta, db, out_tsv, threads=4, evalue=1e-5, max_targets=5, block=0.5):
    """DIAMOND blastx of the reads against the whole family (same settings as the census).

    The read-sequence columns are added when this DIAMOND supports them; older versions
    give the seven columns of version 0.1.0 and the stop, GC and length diagnostics are blank.
    """
    exe = which("diamond", "READLOSS_DIAMOND")
    fields = DIAMOND_FIELDS if seq_fields_supported(exe) else DIAMOND_FIELDS[:7]
    cmd = [exe, "blastx", "-q", reads_fasta, "-d", db, "-e", str(evalue), "-k", str(max_targets),
           "--max-hsps", "1", "-p", str(threads), "-b", str(block), "-c", "1", "--quiet",
           "-f", "6", *fields, "-o", out_tsv]
    subprocess.run(cmd, check=True)
    return out_tsv


# --- nucleotide mode (version 0.3.0) ---
#
# Homeologs of an allopolyploid and recent duplicates can encode the same protein and differ
# only at a few nucleotides, so a protein search cannot tell them apart. In nucleotide mode the
# family is a FASTA of coding sequences, reads are aligned with bwa-mem2 mem -a (all
# alignments), and each alignment becomes a row of the usual hit table:
#   identity  = 100 x (aligned bases - NM) / aligned bases
#   length    = aligned bases / 3, so that --min-aln keeps its meaning in codons
#   bitscore  = AS, the bwa alignment score (one mismatch costs 5 points with default scoring)
#   sstart/send in codons of the coding sequence; the read sequence for GC and length.

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def sam_to_hits(lines):
    """Rows of the hit table from SAM lines (header and unmapped records skipped)."""
    for line in lines:
        if not line or line[0] == "@":
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11 or int(f[1]) & 4:
            continue
        ops = _CIGAR.findall(f[5])
        aligned = sum(int(n) for n, o in ops if o in "M=XI")
        span = sum(int(n) for n, o in ops if o in "M=XDN")
        if not aligned:
            continue
        nm = asc = None
        for t in f[11:]:
            if t.startswith("NM:i:"):
                nm = int(t[5:])
            elif t.startswith("AS:i:"):
                asc = int(t[5:])
        if asc is None:
            continue
        nm = nm or 0
        pos = int(f[3])
        yield (f[0], f[2], round(100 * (aligned - nm) / aligned, 2), aligned // 3, asc,
               (pos - 1) // 3 + 1, (pos - 2 + span) // 3 + 1, "-", f[9])


def write_hits(rows, out_tsv):
    with open(out_tsv, "w") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    return out_tsv


def makedb_nt(cds_fna, out_db):
    """Index a FASTA of coding sequences for the nucleotide mode; returns the indexed FASTA."""
    import shutil as _sh
    exe = which("bwa-mem2", "READLOSS_BWA")
    dst = out_db if out_db.endswith((".fna", ".fa", ".fasta")) else out_db + ".fna"
    if os.path.abspath(dst) != os.path.abspath(cds_fna):
        _sh.copyfile(cds_fna, dst)
    subprocess.run([exe, "index", dst], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return dst


def search_nt(reads, db_fna, out_tsv, threads=4):
    """bwa-mem2 mem -a of the reads (FASTA or FASTQ, single-end) against the coding sequences."""
    exe = which("bwa-mem2", "READLOSS_BWA")
    p = subprocess.Popen([exe, "mem", "-a", "-t", str(threads), db_fna, reads], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, text=True)
    write_hits(sam_to_hits(p.stdout), out_tsv)
    if p.wait() != 0:
        raise subprocess.CalledProcessError(p.returncode, "bwa-mem2 mem")
    return out_tsv


def sam_file_to_hits(path, out_tsv):
    """SAM (text) or BAM/CRAM (through samtools view) of alignments to the coding sequences."""
    if path.endswith((".bam", ".cram")):
        exe = which("samtools", "READLOSS_SAMTOOLS")
        p = subprocess.Popen([exe, "view", path], stdout=subprocess.PIPE, text=True)
        write_hits(sam_to_hits(p.stdout), out_tsv)
        if p.wait() != 0:
            raise subprocess.CalledProcessError(p.returncode, "samtools view")
    else:
        with open(path) as fh:
            write_hits(sam_to_hits(fh), out_tsv)
    return out_tsv


def search_proteins(query_faa, db, out_tsv, threads=4, evalue=1e-5, max_targets=5, block=0.5):
    """DIAMOND blastp of protein fragments against the family, with the settings of ``search``.

    Used by ``readloss selftest``: a fragment of 50 aa stands for a read of 150 nt.
    """
    exe = which("diamond", "READLOSS_DIAMOND")
    cmd = [exe, "blastp", "-q", query_faa, "-d", db, "-e", str(evalue), "-k", str(max_targets),
           "--max-hsps", "1", "-p", str(threads), "-b", str(block), "-c", "1", "--quiet",
           "-f", "6", *DIAMOND_FIELDS[:7], "-o", out_tsv]        # read sequence fields need a translated search
    subprocess.run(cmd, check=True)
    return out_tsv
