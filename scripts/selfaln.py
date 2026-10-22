#!/usr/bin/env python3
"""
selfaln.py - self-alignment submodule for run_pipeline.py.

Called as: selfaln.py [-t THREADS] [--tmpdir DIR] [--fastga-bin PATH] -o out.1aln seqfile

Runs FastGA (https://github.com/thegenemyers/FASTGA) to produce an all-vs-all
self-alignment of the input genome. The output format follows the extension of
the output file name:

    out.1aln (default)  FastGA -v -T<threads> -P<tmpdir> -1:<out> <seqfile>
    out.paf             FastGA -v -T<threads> -P<tmpdir> -paf <seqfile> > <out>
    out.paf.gz          FastGA -v -T<threads> -P<tmpdir> -paf <seqfile> | gzip -c > <out>

The .1aln output is recommended: it is compact and carries the GDB skeleton of the
genome that hapscotch, hapcure, hapcount and alnsplit read directly.
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
    p.add_argument("-o", "--out", required=True,
                   help="output self-alignment path; .1aln (recommended), .paf or .paf.gz")
    p.add_argument("-t", "--threads", type=int, default=8, help="threads [8]")
    p.add_argument("--tmpdir", default=".", help="scratch directory for FastGA -P [.]")
    p.add_argument("--fastga-bin", default="FastGA",
                   help="path to the FastGA executable [FastGA on PATH]")
    return p


def output_format(out: str) -> str:
    """Output format from the output file name: '1aln', 'paf' or 'paf.gz'."""
    if out.endswith(".paf.gz"):
        return "paf.gz"
    if out.endswith(".paf"):
        return "paf"
    if out.endswith(".1aln"):
        return "1aln"
    raise ValueError(f"cannot tell output format from file name: {out} "
                     f"(expect .1aln, .paf or .paf.gz)")


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
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    try:
        output_format(args.out)
    except ValueError as e:
        parser.error(str(e))

    _add_bin_dir_to_path(args.fastga_bin)

    if shutil.which(args.fastga_bin) is None:
        print(f"[selfaln] error: FastGA executable not found: {args.fastga_bin} "
              f"(pass --fastga-bin PATH or add it to PATH)", file=sys.stderr)
        return 1

    base = "{fastga} -v -T{t} -P{tmp}".format(
        fastga=shlex.quote(args.fastga_bin),
        t=args.threads,
        tmp=shlex.quote(args.tmpdir),
    )
    seq = shlex.quote(args.seqfile)
    out = shlex.quote(args.out)
    fmt = output_format(args.out)
    if fmt == "1aln":
        # binary ONEcode output must be written to a file, not streamed
        cmd = f"{base} -1:{out} {seq}"
    elif fmt == "paf.gz":
        # pipefail so that a FastGA failure is not masked by gzip
        cmd = f"set -o pipefail; {base} -paf {seq} | gzip -c > {out}"
    else:
        cmd = f"{base} -paf {seq} > {out}"
    print(f"[selfaln] running: {cmd}", file=sys.stderr)
    return subprocess.run(cmd, shell=True, executable="/bin/bash").returncode


if __name__ == "__main__":
    sys.exit(main())
