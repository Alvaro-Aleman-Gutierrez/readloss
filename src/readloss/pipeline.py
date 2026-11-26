"""External steps: partial download from the SRA and DIAMOND search."""
import os
import re
import shutil
import subprocess
import time

DIAMOND_FIELDS = ["qseqid", "sseqid", "pident", "length", "bitscore", "sstart", "send"]


def which(tool, env=None):
    path = os.environ.get(env) if env else None
    path = path or shutil.which(tool)
    if not path:
        raise FileNotFoundError(f"{tool} not found in PATH" + (f" (or set {env})" if env else ""))
    return path


def fetch(run, spots, out_fasta, retries=3, wait=30, log=None, allow_short=False):
    """Download the first ``spots`` spots of a run as FASTA, mates split.

    The download is accepted only if fastq-dump reports exactly ``spots`` written,
    so a truncated transfer is never analysed. With ``allow_short`` a run that is
    smaller than requested is accepted as it is (check the log, because a broken
    transfer that exits cleanly would look the same).
    """
    exe = which("fastq-dump", "READLOSS_FASTQ_DUMP")
    cmd = [exe, "--fasta", "0", "--split-spot", "--skip-technical", "-X", str(spots), "--stdout", run]
    for i in range(1, retries + 1):
        with open(out_fasta, "w") as out:
            p = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
        if log:
            with open(log, "a") as fh:
                fh.write(p.stderr)
        m = re.search(r"^Written (\d+) spots", p.stderr, re.M)
        if p.returncode == 0 and m and int(m.group(1)) == spots:
            return int(m.group(1))
        if allow_short and m and p.returncode == 0 and 0 < int(m.group(1)) < spots:
            return int(m.group(1))
        os.remove(out_fasta)
        if i < retries:
            time.sleep(wait)
    raise RuntimeError(f"download of {run} failed after {retries} attempts")


def fetch_search_stream(run, spots, db, out_tsv, threads=4, retries=3, wait=30, log=None,
                        allow_short=False):
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
        dump = subprocess.Popen(
            ["sh", "-c", f'exec "$0" --fasta 0 --split-spot --skip-technical -X {spots} '
                         f'--stdout {run} > "$1"', fq, fifo],
            stderr=subprocess.PIPE, text=True)
        try:
            search(fifo, db, out_tsv, threads=threads)
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
        m = re.search(r"^Written (\d+) spots", err, re.M)
        n = int(m.group(1)) if m else 0
        if ok_search and dump.returncode == 0 and (n == spots or (allow_short and 0 < n < spots)):
            return n
        if os.path.exists(out_tsv):
            os.remove(out_tsv)
        if i < retries:
            time.sleep(wait)
    raise RuntimeError(f"streamed download and search of {run} failed after {retries} attempts")


def makedb(family_faa, out_db):
    exe = which("diamond", "READLOSS_DIAMOND")
    subprocess.run([exe, "makedb", "--in", family_faa, "-d", out_db, "--quiet"], check=True)
    return out_db if out_db.endswith(".dmnd") else out_db + ".dmnd"


def search(reads_fasta, db, out_tsv, threads=4, evalue=1e-5, max_targets=5, block=0.5):
    """DIAMOND blastx of the reads against the whole family (same settings as the census)."""
    exe = which("diamond", "READLOSS_DIAMOND")
    cmd = [exe, "blastx", "-q", reads_fasta, "-d", db, "-e", str(evalue), "-k", str(max_targets),
           "--max-hsps", "1", "-p", str(threads), "-b", str(block), "-c", "1", "--quiet",
           "-f", "6", *DIAMOND_FIELDS, "-o", out_tsv]
    subprocess.run(cmd, check=True)
    return out_tsv


def search_proteins(query_faa, db, out_tsv, threads=4, evalue=1e-5, max_targets=5, block=0.5):
    """DIAMOND blastp of protein fragments against the family, with the settings of ``search``.

    Used by ``readloss selftest``: a fragment of 50 aa stands for a read of 150 nt.
    """
    exe = which("diamond", "READLOSS_DIAMOND")
    cmd = [exe, "blastp", "-q", query_faa, "-d", db, "-e", str(evalue), "-k", str(max_targets),
           "--max-hsps", "1", "-p", str(threads), "-b", str(block), "-c", "1", "--quiet",
           "-f", "6", *DIAMOND_FIELDS, "-o", out_tsv]
    subprocess.run(cmd, check=True)
    return out_tsv
