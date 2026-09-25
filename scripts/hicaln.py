#!/usr/bin/env python3
"""
hicaln.py - HiC read alignment submodule for run_pipeline.py.

Called as: hicaln.py [-t THREADS] -o out_prefix seqfile hicfile

Aligns one HiC input against seqfile using minibwa (--hic mode) or bwa-mem2
and produces a name-sorted, duplicate-marked BAM at <out_prefix>.bam. All working
files (reference index, intermediates) are written alongside <out_prefix> - the
caller is expected to point out_prefix at a dedicated directory (e.g.
data/hicalns/hicaln.0) shared across all HiC inputs, since the reference index
is built once in that directory and reused (skipped if already present) for
every subsequent call.

hicfile may be:
  - a single .bam/.cram file (reads are re-extracted to FASTA via `samtools
    fasta` before alignment)
  - a single fasta/fastq[.gz] file
  - two comma-separated mate fasta/fastq[.gz] files, e.g. R1.fq.gz,R2.fq.gz

Pipeline per input, all under the out_prefix's directory:
  [once]  aligner index <seqfile> <dir>/ref
  bam/cram:  samtools fasta -F0xB00 -n <in> | aligner map/mem ref - |
             samtools fixmate -mpu - - | samtools sort --write-index -l1
             -T <out_prefix>.tmp -o <out_prefix>.srt.bam -
  fasta/q:   aligner map/mem ref <in1> [in2] | samtools fixmate -mpu
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


def _index_genome(seqfile: Path, workdir: Path, threads: int, aligner: str, aligner_bin: str) -> Path:
    ref_prefix = workdir / "ref"
    exts = [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"] if aligner == "bwamem2" else [".l2b", ".mbw"]
    idx_files = [Path(f"{ref_prefix}{ext}") for ext in exts]
    if all(f.exists() for f in idx_files):
        _log("index", f"skip (already indexed: {ref_prefix})")
        return ref_prefix
    if aligner == "bwamem2":
        _run_stage("index", [aligner_bin, "index", "-p", str(ref_prefix), str(seqfile)])
    else:
        index_threads = min(threads, 4)
        _run_stage("index", [aligner_bin, "index", f"-t{index_threads}", str(seqfile), str(ref_prefix)])
    if not all(f.exists() for f in idx_files):
        raise HicAlnError(f"[index] completed but expected output(s) missing: {ref_prefix}")
    return ref_prefix


def _align(hicfile: str, ref_prefix: Path, out_prefix: str, threads: int,
           aligner: str, aligner_bin: str, samtools_bin: str) -> str:
    srt_bam = f"{out_prefix}.srt.bam"
    tmp_prefix = f"{out_prefix}.tmp"
    mode = "mem -5SP" if aligner == "bwamem2" else "map --hic"
    map_cmd = f"{shlex.quote(aligner_bin)} {mode} -t {threads} {shlex.quote(str(ref_prefix))}"
    if _is_bam_cram(hicfile):
        if aligner == "bwamem2":
            map_cmd += " -p"
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
        if aligner == "bwamem2" and len(inputs) == 1:
            map_cmd += " -p"
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
                        "(including the shared reference index) are written next to it")
    p.add_argument("-a", "--aligner", choices=["bwamem2", "minibwa"], default="bwamem2",
                   help="aligner program to use [bwamem2]")
    p.add_argument("-t", "--threads", type=int, default=8, help="threads [8]")
    p.add_argument("--bwamem2-bin", default="bwa-mem2", help="path to the bwa-mem2 executable [bwa-mem2]")
    p.add_argument("--minibwa-bin", default="minibwa", help="path to the minibwa executable [minibwa]")
    p.add_argument("--samtools-bin", default="samtools", help="path to the samtools executable [samtools]")
    return p


def main(argv=None) -> int:
    args = build_arg_parser().parse_args(argv)

    aligner_bin = args.bwamem2_bin if args.aligner == "bwamem2" else args.minibwa_bin
    missing = [t for t in (aligner_bin, args.samtools_bin) if shutil.which(t) is None]
    if missing:
        _log("preflight", f"error: required tool(s) not found on PATH: {', '.join(missing)}")
        return 1

    workdir = Path(args.out_prefix).parent
    workdir.mkdir(parents=True, exist_ok=True)

    try:
        ref_prefix = _index_genome(Path(args.seqfile), workdir, args.threads, args.aligner, aligner_bin)
        srt_bam = _align(args.hicfile, ref_prefix, args.out_prefix, args.threads,
                          args.aligner, aligner_bin, args.samtools_bin)
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
