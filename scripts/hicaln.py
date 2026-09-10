#!/usr/bin/env python3
"""
hicaln.py - HiC read alignment submodule for run_pipeline.py.

Called as: hicaln.py [-t THREADS] -o out_prefix seqfile hicfile

Aligns one HiC input against seqfile using minibwa (--hic mode) and produces a
name-sorted, duplicate-marked BAM at <out_prefix>.bam. All working files
(reference index, intermediates) are written alongside <out_prefix> - the
caller is expected to point out_prefix at a dedicated directory (e.g.
data/hicalns/hicaln.0) shared across all HiC inputs, since the minibwa
reference index (<dir>/ref.l2b + <dir>/ref.mbw) is built once in that
directory and reused (skipped if already present) for every subsequent call.

hicfile may be:
  - a single .bam/.cram file (reads are re-extracted to FASTA via `samtools
    fasta` before alignment)
  - a single fasta/fastq[.gz] file
  - two comma-separated mate fasta/fastq[.gz] files, e.g. R1.fq.gz,R2.fq.gz

Pipeline per input, all under the out_prefix's directory:
  [once]  minibwa index -t<=4> <seqfile> <dir>/ref
  bam/cram:  samtools fasta -F0xB00 -n <in> | minibwa map --hic -t<T> ref - |
             samtools fixmate -mpu - - | samtools sort --write-index -l1
             -T <out_prefix>.tmp -o <out_prefix>.srt.bam -
  fasta/q:   minibwa map --hic -t<T> ref <in1> [in2] | samtools fixmate -mpu
             - - | samtools sort --write-index -l1 -T <out_prefix>.tmp
             -o <out_prefix>.srt.bam -
  samtools markdup --write-index -c -@<T> -T <out_prefix>.mkdup.tmp
      <out_prefix>.srt.bam <out_prefix>.mkdup.bam
  samtools sort -N -@<T> -T <out_prefix>.tmp -o <out_prefix>.bam
      <out_prefix>.mkdup.bam
  (then .srt.bam/.mkdup.bam and their indexes are removed)
"""
import argparse
import glob
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


class HicAlnError(RuntimeError):
    """Raised when a stage of the HiC alignment pipeline fails."""


def _log(stage: str, msg: str) -> None:
    print(f"[hicaln] [{stage}] {msg}", file=sys.stderr, flush=True)


def _run_stage(stage: str, cmd, *, shell: bool = False) -> None:
    printable = cmd if isinstance(cmd, str) else " ".join(shlex.quote(str(c)) for c in cmd)
    _log(stage, f"running: {printable}")
    try:
        proc = subprocess.run(cmd, shell=shell)
    except OSError as e:
        raise HicAlnError(f"[{stage}] could not launch command: {printable}\n{e}") from e
    if proc.returncode != 0:
        raise HicAlnError(f"[{stage}] command failed (exit {proc.returncode}): {printable}")


def _is_bam_cram(path: str) -> bool:
    return Path(path).suffix.lower() in (".bam", ".cram")


def _index_genome(seqfile: Path, workdir: Path, threads: int, minibwa_bin: str) -> Path:
    ref_prefix = workdir / "ref"
    l2b = Path(f"{ref_prefix}.l2b")
    mbw = Path(f"{ref_prefix}.mbw")
    if l2b.exists() and mbw.exists():
        _log("index", f"skip (already indexed: {l2b}, {mbw})")
        return ref_prefix
    index_threads = min(threads, 4)  # indexing itself doesn't benefit past ~4 threads
    _run_stage("index", [minibwa_bin, "index", f"-t{index_threads}", str(seqfile), str(ref_prefix)])
    if not (l2b.exists() and mbw.exists()):
        raise HicAlnError(f"[index] completed but expected output(s) missing: {l2b}, {mbw}")
    return ref_prefix


def _align(hicfile: str, ref_prefix: Path, out_prefix: str, threads: int,
           minibwa_bin: str, samtools_bin: str) -> str:
    srt_bam = f"{out_prefix}.srt.bam"
    tmp_prefix = f"{out_prefix}.tmp"
    map_cmd = f"{shlex.quote(minibwa_bin)} map --hic -t {threads} {shlex.quote(str(ref_prefix))}"
    if _is_bam_cram(hicfile):
        cmd = (
            f"{shlex.quote(samtools_bin)} fasta -F0xB00 -n {shlex.quote(hicfile)} | "
            f"{map_cmd} - | "
            f"{shlex.quote(samtools_bin)} fixmate -mpu - - | "
            f"{shlex.quote(samtools_bin)} sort --write-index -l1 -T {shlex.quote(tmp_prefix)} "
            f"-o {shlex.quote(srt_bam)} -"
        )
    else:
        inputs = [s.strip() for s in hicfile.split(",")]
        if len(inputs) > 2:
            raise HicAlnError(f"[align] expected at most 2 comma-separated mate files, "
                               f"got {len(inputs)}: {hicfile}")
        quoted_inputs = " ".join(shlex.quote(i) for i in inputs)
        cmd = (
            f"{map_cmd} {quoted_inputs} | "
            f"{shlex.quote(samtools_bin)} fixmate -mpu - - | "
            f"{shlex.quote(samtools_bin)} sort --write-index -l1 -T {shlex.quote(tmp_prefix)} "
            f"-o {shlex.quote(srt_bam)} -"
        )
    _run_stage("align", cmd, shell=True)
    if not Path(srt_bam).exists():
        raise HicAlnError(f"[align] completed but expected output missing: {srt_bam}")
    return srt_bam


def _markdup(srt_bam: str, out_prefix: str, threads: int, samtools_bin: str) -> str:
    mkdup_bam = f"{out_prefix}.mkdup.bam"
    mkdup_tmp = f"{out_prefix}.mkdup.tmp"
    _run_stage("markdup", [samtools_bin, "markdup", "--write-index", "-c", "-@", str(threads),
                            "-T", mkdup_tmp, srt_bam, mkdup_bam])
    if not Path(mkdup_bam).exists():
        raise HicAlnError(f"[markdup] completed but expected output missing: {mkdup_bam}")
    return mkdup_bam


def _name_sort(mkdup_bam: str, out_prefix: str, threads: int, samtools_bin: str) -> str:
    out_bam = f"{out_prefix}.bam"
    tmp_prefix = f"{out_prefix}.tmp"
    _run_stage("sort", [samtools_bin, "sort", "-N", "-@", str(threads),
                         "-T", tmp_prefix, "-o", out_bam, mkdup_bam])
    if not Path(out_bam).exists():
        raise HicAlnError(f"[sort] completed but expected output missing: {out_bam}")
    return out_bam


def _cleanup(*prefixed_paths: str) -> None:
    for p in prefixed_paths:
        for f in glob.glob(f"{p}*"):
            _log("cleanup", f"removing {f}")
            try:
                Path(f).unlink()
            except OSError as e:
                _log("cleanup", f"warning: could not remove {f}: {e}")


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("seqfile", help="genome sequence file (fasta[.gz]) to index and align against")
    p.add_argument("hicfile", help="raw HiC input: a single BAM/CRAM file, a single "
                                   "fasta/fastq[.gz] file, or two comma-separated mate "
                                   "fasta/fastq[.gz] files")
    p.add_argument("-o", "--out-prefix", required=True, metavar="PREFIX",
                   help="output prefix; final alignment is PREFIX.bam, all working files "
                        "(including the shared minibwa index) are written next to it")
    p.add_argument("-t", "--threads", type=int, default=8, help="threads [8]")
    p.add_argument("--minibwa-bin", default="minibwa", help="path to the minibwa executable [minibwa]")
    p.add_argument("--samtools-bin", default="samtools", help="path to the samtools executable [samtools]")
    return p


def main(argv=None) -> int:
    args = build_arg_parser().parse_args(argv)

    missing = [t for t in (args.minibwa_bin, args.samtools_bin) if shutil.which(t) is None]
    if missing:
        _log("preflight", f"error: required tool(s) not found on PATH: {', '.join(missing)}")
        return 1

    workdir = Path(args.out_prefix).parent
    workdir.mkdir(parents=True, exist_ok=True)

    try:
        ref_prefix = _index_genome(Path(args.seqfile), workdir, args.threads, args.minibwa_bin)
        srt_bam = _align(args.hicfile, ref_prefix, args.out_prefix, args.threads,
                          args.minibwa_bin, args.samtools_bin)
        mkdup_bam = _markdup(srt_bam, args.out_prefix, args.threads, args.samtools_bin)
        out_bam = _name_sort(mkdup_bam, args.out_prefix, args.threads, args.samtools_bin)
        _cleanup(srt_bam, mkdup_bam)
    except HicAlnError as e:
        _log("error", str(e))
        return 1

    _log("done", out_bam)
    return 0


if __name__ == "__main__":
    sys.exit(main())
