#!/usr/bin/env python3
"""
selfaln.py - self-alignment submodule for run_pipeline.py.

Called as: selfaln.py [-t THREADS] [--tmpdir DIR] -o out.paf seqfile

Placeholder: no aligner is wired in yet. Edit the TODO block below to invoke
your preferred self-alignment tool (e.g. minimap2 ava-ont/ava-pb) and write a
PAF file to the path given by -o.
"""
import argparse
import subprocess
import sys
from pathlib import Path


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("seqfile", help="genome sequence file (fasta[.gz])")
    p.add_argument("-o", "--out", required=True, help="output self-alignment PAF path")
    p.add_argument("-t", "--threads", type=int, default=8, help="threads [8]")
    p.add_argument("--tmpdir", default=".", help="scratch directory [.]")
    return p


def main(argv=None) -> int:
    args = build_arg_parser().parse_args(argv)

    # TODO: replace with a real self-alignment command, e.g.:
    #   subprocess.run(
    #       f"minimap2 -x ava-ont -t {args.threads} {args.seqfile} {args.seqfile} "
    #       f"> {args.out}", shell=True, check=True)
    print(
        "[selfaln] not implemented: edit selfaln.py's TODO block to run your "
        f"preferred self-aligner, writing PAF output to {args.out}",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
