#!/usr/bin/env python3
"""
hicaln.py - HiC read alignment submodule for run_pipeline.py.

Called as: hicaln.py [-t THREADS] [--tmpdir DIR] -o out.bam seqfile hicfile

Placeholder: no aligner is wired in yet. Edit the TODO block below to invoke
your preferred HiC aligner (e.g. bwa mem + samtools sort -n) and write a
name-sorted BAM (or BED/PA5) file to the path given by -o.
"""
import argparse
import subprocess
import sys
from pathlib import Path


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("seqfile", help="genome sequence file (fasta[.gz])")
    p.add_argument("hicfile", help="raw HiC fastq/a file to align")
    p.add_argument("-o", "--out", required=True, help="output alignment path (e.g. .bam)")
    p.add_argument("-t", "--threads", type=int, default=8, help="threads [8]")
    p.add_argument("--tmpdir", default=".", help="scratch directory [.]")
    return p


def main(argv=None) -> int:
    args = build_arg_parser().parse_args(argv)

    # TODO: replace with a real HiC alignment command, e.g.:
    #   subprocess.run(
    #       f"bwa mem -t {args.threads} {args.seqfile} {args.hicfile} | "
    #       f"samtools sort -n -o {args.out} -", shell=True, check=True)
    print(
        "[hicaln] not implemented: edit hicaln.py's TODO block to run your "
        f"preferred HiC aligner on {args.hicfile}, writing output to {args.out}",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
