#!/usr/bin/env python3
"""
selfaln.py - self-alignment submodule for run_pipeline.py.

Called as: selfaln.py [-t THREADS] [--tmpdir DIR] [--fastga-bin PATH] -o out.paf seqfile

Runs FastGA (https://github.com/thegenemyers/FASTGA) to produce an all-vs-all
self-alignment PAF of the input genome:

    FastGA -v -T<threads> -P<tmpdir> -paf <seqfile> > <out>
"""
import argparse
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("seqfile", help="genome sequence file (fasta[.gz])")
    p.add_argument("-o", "--out", required=True, help="output self-alignment PAF path")
    p.add_argument("-t", "--threads", type=int, default=8, help="threads [8]")
    p.add_argument("--tmpdir", default=".", help="scratch directory for FastGA -P [.]")
    p.add_argument("--fastga-bin", default="FastGA",
                   help="path to the FastGA executable [FastGA on PATH]")
    return p


def _add_bin_dir_to_path(fastga_bin: str) -> None:
    """FastGA internally calls several sibling tools (e.g. GIXmake, ALNtoPAF) - if a
    path (not just a bare name) is given, add its directory to PATH so those are
    found too, since they live alongside the FastGA executable itself."""
    bin_dir = os.path.dirname(fastga_bin)
    if not bin_dir:
        return
    path_dirs = os.environ.get("PATH", "").split(os.pathsep)
    if bin_dir not in path_dirs:
        os.environ["PATH"] = os.pathsep.join([bin_dir] + path_dirs)


def main(argv=None) -> int:
    args = build_arg_parser().parse_args(argv)

    _add_bin_dir_to_path(args.fastga_bin)

    if shutil.which(args.fastga_bin) is None:
        print(f"[selfaln] error: FastGA executable not found: {args.fastga_bin} "
              f"(pass --fastga-bin PATH or add it to PATH)", file=sys.stderr)
        return 1

    cmd = "{fastga} -v -T{t} -P{tmp} -paf {seq} > {out}".format(
        fastga=shlex.quote(args.fastga_bin),
        t=args.threads,
        tmp=shlex.quote(args.tmpdir),
        seq=shlex.quote(args.seqfile),
        out=shlex.quote(args.out),
    )
    print(f"[selfaln] running: {cmd}", file=sys.stderr)
    return subprocess.run(cmd, shell=True).returncode


if __name__ == "__main__":
    sys.exit(main())
