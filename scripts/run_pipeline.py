#!/usr/bin/env python3
"""
run_pipeline.py - end-to-end driver for the HapScotch haplotyping pipeline.

Stages (each is independently resumable - see Checkpointing below):
  1. seq_index       seqtools idx                   -> genome sequence index (.idx)
  2. self_align      selfaln.py                     -> self-alignment PAF
  3. hic_align       hicaln.py                      -> per-HICFILE alignment (BAM/...)
  4. hic_convert     hictools convert               -> merged HiC BIN
  5. contig_ec       hapcure                        -> contig error-correction AGP
  6. hapscotch       hapscotch                      -> haplotype phased scaffolds
  7. yahs_scaffold   seqtools/yahs/hictools/hicmap  -> HiC-scaffolded haplotypes
  8. collect_results copy final files               -> <outdir>/4.results

Checkpointing: every step declares the output file(s) it must produce. On start-up
the pipeline checks, for each step in order, whether its declared outputs already
exist; if so the step is skipped (logged as "[skip]"). Use --resume to continue
a previous run from wherever it stopped (completed steps are skipped, nothing is
rerun unnecessarily); --force to ignore existing outputs and rerun every step from
scratch; or --force-from STEP to rerun a step and everything after it. If OUTDIR
already exists, one of --resume/--force/--force-from is required (otherwise the
pipeline refuses to run); conversely, those options require OUTDIR to already
exist. The parameters of the run that created OUTDIR are recorded in
OUTDIR/data/CMD; --resume/--force/--force-from are refused if the current
parameters don't match that record.
"""
import argparse
import glob
import json
import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

# --------------------------------------------------------------------------- #
# errors / logging
# --------------------------------------------------------------------------- #

class PipelineError(RuntimeError):
    """Raised when a step fails or its expected outputs are missing after running."""


def _log(msg: str) -> None:
    print(f"[run_pipeline] {msg}", file=sys.stderr, flush=True)


# sibling python submodules - always resolved next to this script, not user-configurable
SCRIPT_DIR = Path(__file__).resolve().parent
HICMAP_PY = SCRIPT_DIR / "hicmap.py"
SELFALN_PY = SCRIPT_DIR / "selfaln.py"
HICALN_PY = SCRIPT_DIR / "hicaln.py"

# compiled binaries (seqtools/hictools/hapscotch/hapcure) live one level up
REPO_ROOT = SCRIPT_DIR.parent


def _add_repo_root_to_path() -> None:
    """so bare tool names (seqtools, hictools, hapscotch, hapcure) resolve
    without needing --*-bin, even when run_pipeline.py is invoked from elsewhere"""
    root = str(REPO_ROOT)
    path_dirs = os.environ.get("PATH", "").split(os.pathsep)
    if root not in path_dirs:
        os.environ["PATH"] = os.pathsep.join([root] + path_dirs)


# --------------------------------------------------------------------------- #
# shared pipeline context (paths, binaries, resources)
# --------------------------------------------------------------------------- #

@dataclass
class Ctx:
    seqfile: Path
    outdir: Path
    threads: int
    verbose: int

    seqtools: str = "seqtools"
    hictools: str = "hictools"
    hapscotch_bin: str = "hapscotch"
    hapcure_bin: str = "hapcure"
    yahs_bin: str = "yahs"

    @property
    def datadir(self) -> Path:
        return self.outdir / "data"

    @property
    def logdir(self) -> Path:
        return self.outdir / "logs"

    @property
    def idxfile(self) -> Path:
        return self.datadir / "seq.idx"

    @property
    def hapcure_dir(self) -> Path:
        return self.outdir / "1.hapcure"

    @property
    def hapscotch_dir(self) -> Path:
        return self.outdir / "2.hapscotch"

    @property
    def yahs_dir(self) -> Path:
        return self.outdir / "3.yahs"

    @property
    def results_dir(self) -> Path:
        return self.outdir / "4.results"

    def ensure_dirs(self) -> None:
        for d in (self.outdir, self.datadir, self.logdir):
            d.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------- #
# command execution
# --------------------------------------------------------------------------- #

def run_cmd(cmd, log_path: Path, *, shell: bool = False, cwd: Optional[Path] = None) -> None:
    """Run `cmd` (list of args, or a string if shell=True), tee-ing combined
    stdout/stderr to log_path. Raises PipelineError on non-zero exit or if the
    executable itself cannot be found/launched."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    printable = cmd if isinstance(cmd, str) else " ".join(shlex.quote(str(c)) for c in cmd)
    _log(f"running: {printable}")
    _log(f"  log: {log_path}")
    with open(log_path, "w") as lf:
        lf.write(f"# CMD: {printable}\n")
        lf.flush()
        try:
            proc = subprocess.run(
                cmd, shell=shell, cwd=cwd, stdout=lf, stderr=subprocess.STDOUT
            )
        except OSError as e:
            raise PipelineError(f"could not launch command: {printable}\n{e}") from e
    if proc.returncode != 0:
        tail = _tail(log_path, 30)
        raise PipelineError(
            f"command failed (exit {proc.returncode}): {printable}\n"
            f"--- last {len(tail)} lines of {log_path} ---\n" + "".join(tail)
        )


def _tail(path: Path, n: int) -> list:
    try:
        with open(path) as f:
            lines = f.readlines()
        return lines[-n:]
    except OSError:
        return []


def require_outputs(outputs) -> None:
    missing = [str(o) for o in outputs if not Path(o).exists()]
    if missing:
        raise PipelineError(
            "step completed but expected output(s) missing: " + ", ".join(missing)
        )


# --------------------------------------------------------------------------- #
# preflight - fail fast if a tool a step needs isn't available, before running anything
# --------------------------------------------------------------------------- #

def _check_binary(name: str, target: str) -> Optional[str]:
    if shutil.which(target) is None:
        return f"{name}: '{target}' not found on PATH (or not executable)"
    return None


def _check_script(name: str, target: Path) -> Optional[str]:
    if not target.is_file():
        return f"{name}: expected script not found: {target}"
    return None


def preflight_check(ctx: Ctx, steps: list) -> None:
    """Verify every external tool/script needed by the steps that are actually
    applicable to this run (Step.applicable()) can be found, before anything
    runs - e.g. yahs/hictools-prepare/hicmap.py are only required when the
    yahs_scaffold step is actually going to run."""
    checks_by_step = {
        "seq_index":       [("seqtools", ctx.seqtools)],
        "self_align":      [("selfaln.py", SELFALN_PY)],
        "hic_align":       [("hicaln.py", HICALN_PY)],
        "hic_convert":     [("hictools", ctx.hictools)],
        "contig_ec":       [("hapcure", ctx.hapcure_bin)],
        "hapscotch":       [("hapscotch", ctx.hapscotch_bin)],
        "yahs_scaffold":   [("seqtools", ctx.seqtools), ("yahs", ctx.yahs_bin),
                            ("hictools", ctx.hictools), ("hicmap.py", HICMAP_PY)],
        "collect_results": [("seqtools", ctx.seqtools)],
    }

    errors = []
    seen = set()
    for step in steps:
        if not step.applicable():
            continue
        for name, target in checks_by_step.get(step.name, []):
            key = (name, str(target))
            if key in seen:
                continue
            seen.add(key)
            err = _check_script(name, target) if isinstance(target, Path) \
                else _check_binary(name, target)
            if err:
                errors.append(f"[{step.name}] {err}")

    if errors:
        raise PipelineError(
            "required tool(s) not found - aborting before running any step:\n  "
            + "\n  ".join(errors)
        )


# --------------------------------------------------------------------------- #
# step framework (output-file-existence based checkpointing)
# --------------------------------------------------------------------------- #

@dataclass
class Step:
    name: str
    outputs: Callable[[], list]   # -> list of Path that must exist when done
    run: Callable[[], None]
    applicable: Callable[[], bool] = field(default=lambda: True)  # step even needed?


def run_step(step: Step, *, force: bool, force_from: bool) -> None:
    if not step.applicable():
        _log(f"[skip] {step.name} (not applicable)")
        return
    outs = step.outputs()
    if not force and not force_from and outs and all(Path(o).exists() for o in outs):
        _log(f"[skip] {step.name} (outputs already present)")
        return
    _log(f"[run ] {step.name}")
    step.run()
    if outs:
        require_outputs(outs)
    _log(f"[done] {step.name}")


# --------------------------------------------------------------------------- #
# step implementations
# --------------------------------------------------------------------------- #

def step_seq_index(ctx: Ctx) -> Step:
    def _run():
        run_cmd(
            [ctx.seqtools, "idx", "-o", str(ctx.idxfile), str(ctx.seqfile)],
            ctx.logdir / "seq_index.log",
        )
    return Step("seq_index", lambda: [ctx.idxfile], _run)


def step_self_align(ctx: Ctx, args) -> Step:
    out = ctx.datadir / "seqaln.paf"

    def _resolved_path() -> Path:
        return Path(args.seq_aln) if args.seq_aln else out

    def _applicable() -> bool:
        return args.seq_aln is None

    def _run():
        run_cmd(
            [sys.executable, str(SELFALN_PY), str(ctx.seqfile), "-o", str(out),
             "-t", str(ctx.threads), "--tmpdir", str(ctx.datadir),
             "--fastga-bin", args.fastga_bin],
            ctx.logdir / "self_align.log",
        )

    return Step("self_align", lambda: [out], _run, applicable=_applicable), _resolved_path


def step_hic_align(ctx: Ctx, args) -> tuple:
    """Aligns each --hic-file input; returns (Step, resolved_hicalns_fn)."""
    hicfiles = args.hic_file or []
    generated = [ctx.datadir / f"hicaln.{i}.bam" for i in range(len(hicfiles))]

    def _applicable() -> bool:
        return len(hicfiles) > 0

    def _run():
        for src, out in zip(hicfiles, generated):
            run_cmd(
                [sys.executable, str(HICALN_PY), str(ctx.seqfile), str(src), "-o", str(out),
                 "-t", str(ctx.threads), "--tmpdir", str(ctx.datadir)],
                ctx.logdir / f"hic_align.{out.stem}.log",
            )

    def _resolved_hicalns() -> list:
        # user-supplied pre-aligned files + freshly generated ones
        return list(args.hic_aln or []) + generated

    return Step("hic_align", lambda: list(generated), _run, applicable=_applicable), _resolved_hicalns


def step_hic_convert(ctx: Ctx, resolved_hicalns_fn) -> tuple:
    out = ctx.datadir / "hicaln.bin"

    def _is_passthrough() -> bool:
        # a single already-BIN input is passed through unchanged, no new file written
        hicalns = resolved_hicalns_fn()
        return len(hicalns) == 1 and str(hicalns[0]).endswith(".bin")

    def _applicable() -> bool:
        return len(resolved_hicalns_fn()) > 0

    def _run():
        hicalns = resolved_hicalns_fn()
        run_cmd(
            [ctx.hictools, "convert", "-o", str(ctx.datadir / "hicaln")]
            + [str(ctx.idxfile)] + [str(h) for h in hicalns],
            ctx.logdir / "hic_convert.log",
        )

    def _outputs() -> list:
        if not resolved_hicalns_fn() or _is_passthrough():
            return []
        return [out]

    def _resolved_hicbin() -> Optional[Path]:
        hicalns = resolved_hicalns_fn()
        if not hicalns:
            return None
        if out.exists():
            return out
        # hictools convert passed a single already-BIN input straight through
        if _is_passthrough():
            return Path(hicalns[0])
        return None

    def _will_have_hic() -> bool:
        # decidable purely from CLI inputs - safe to call before hic_convert has run,
        # unlike _resolved_hicbin() which depends on files that may not exist yet
        return len(resolved_hicalns_fn()) > 0

    return Step("hic_convert", _outputs, _run,
                applicable=_applicable), _resolved_hicbin, _will_have_hic


def step_contig_ec(ctx: Ctx, args, resolved_hicbin_fn, will_have_hic_fn) -> tuple:
    out = ctx.hapcure_dir / "ctg.ec.agp"

    def _applicable() -> bool:
        return args.contig_ec and will_have_hic_fn()

    def _run():
        ctx.hapcure_dir.mkdir(parents=True, exist_ok=True)
        run_cmd(
            [ctx.hapcure_bin, "-o", str(ctx.hapcure_dir / "ctg"),
             str(ctx.idxfile), str(resolved_hicbin_fn())]
            + args.hapcure_opt,
            ctx.logdir / "contig_ec.log",
        )

    def _resolved_agpec() -> Optional[Path]:
        return out if _applicable() else None

    return Step("contig_ec", lambda: [out], _run, applicable=_applicable), _resolved_agpec


def step_hapscotch(ctx: Ctx, args, seqaln_fn, resolved_hicbin_fn, resolved_agpec_fn,
                    run_yahs_fn) -> Step:
    prefix = ctx.hapscotch_dir / "haps"
    out = Path(f"{prefix}.scfs.agp")

    def _run():
        ctx.hapscotch_dir.mkdir(parents=True, exist_ok=True)
        opts = []
        if args.busco:
            opts += ["-g", str(args.busco)]
        if args.ploidy is not None:
            opts += ["-p", str(args.ploidy)]
        if args.max_ploidy is not None:
            opts += ["--max-ploidy", str(args.max_ploidy)]
        if args.min_ext is not None:
            opts += ["-B", str(args.min_ext)]
        if args.min_qual is not None:
            opts += ["-q", str(args.min_qual)]
        opts += ["-t", str(ctx.threads)]
        agpec = resolved_agpec_fn()
        if agpec is not None:
            opts += ["-a", str(agpec)]
        hicbin = resolved_hicbin_fn()
        if hicbin is not None:
            opts += ["-c", str(hicbin), "--file-type", "BIN"]
        if run_yahs_fn():
            opts += ["-Y"]
        opts += args.hapscotch_opt
        run_cmd(
            [ctx.hapscotch_bin] + opts + ["-o", str(prefix), str(ctx.idxfile), str(seqaln_fn())],
            ctx.logdir / "hapscotch.log",
        )

    return Step("hapscotch", lambda: [out], _run)


def step_yahs_scaffold(ctx: Ctx, args, resolved_hicbin_fn, run_yahs_fn) -> Step:
    hap_prefix = ctx.hapscotch_dir / "haps"
    final_agp = ctx.yahs_dir / "allhaps.scf.agp"

    def _applicable() -> bool:
        return run_yahs_fn()

    def _run():
        ctx.yahs_dir.mkdir(parents=True, exist_ok=True)
        bbseq_agp = Path(f"{hap_prefix}.bbseq.agp")
        bbscf_agp = Path(f"{hap_prefix}.bbscf.agp")
        bbpos_txt = Path(f"{hap_prefix}.bbpos.txt")
        bbscf_hic_bin = Path(f"{hap_prefix}.bbscf-hic.bin")
        seq_hic_bin = resolved_hicbin_fn()
        for needed in (bbseq_agp, bbscf_agp, bbpos_txt, bbscf_hic_bin, seq_hic_bin):
            if not needed.exists():
                raise PipelineError(
                    f"YaHS branch requested but required hapscotch output missing: {needed}\n"
                    "(hapscotch only writes these when run with -c/HiC data and -Y)"
                )

        bbseq_fa = ctx.yahs_dir / "haps.bbseq.fa"
        run_cmd(
            [ctx.seqtools, "seq", "-o", str(bbseq_fa), str(bbseq_agp), str(ctx.seqfile)],
            ctx.logdir / "yahs.seqtools_seq.log",
        )

        bbseq_idx = ctx.yahs_dir / "haps.bbseq.fa.fai" # yahs need .fai not .idx
        run_cmd(
            [ctx.seqtools, "idx", "-o", str(bbseq_idx), str(bbseq_fa)],
            ctx.logdir / "yahs.seqtools_idx.log",
        )

        yahs_prefix = ctx.yahs_dir / "haps.bbseq"
        run_cmd(
            [ctx.yahs_bin, "-a", str(bbscf_agp), "--no-contig-ec", "--no-scaffold-ec",
             "-o", str(yahs_prefix)] + args.yahs_opt + [str(bbseq_fa), str(bbscf_hic_bin)],
            ctx.logdir / "yahs.log",
        )
        yahs_scaffolds = _find_yahs_scaffolds_agp(yahs_prefix)

        run_cmd(
            [ctx.seqtools, "hap", "-o", str(final_agp), str(yahs_scaffolds), str(bbpos_txt)],
            ctx.logdir / "yahs.seqtools_hap.log",
        )

        hic_txt = ctx.yahs_dir / "allhaps.scf.hic.txt"
        run_cmd(
            [ctx.hictools, "prepare", "-a", str(final_agp), "-n", "2000",
             "-o", str(hic_txt)] + args.hictools_prepare_opt
            + [str(seq_hic_bin), str(ctx.idxfile)],
            ctx.logdir / "yahs.hictools_prepare.log",
        )

        hic_png = ctx.yahs_dir / "allhaps.scf.hic.png"
        hic_pdf = ctx.yahs_dir / "allhaps.scf.hic.pdf"
        run_cmd(
            [sys.executable, str(HICMAP_PY), "--png", str(hic_png), "--pdf", str(hic_pdf)]
            + args.hicmap_opt + [str(hic_txt)],
            ctx.logdir / "yahs.hicmap.log",
        )

    return Step("yahs_scaffold", lambda: [final_agp], _run, applicable=_applicable)


def _find_yahs_scaffolds_agp(prefix: Path) -> Path:
    """yahs's exact output-file naming can vary by build; probe the common
    patterns instead of hard-coding one."""
    candidates = [
        Path(f"{prefix}_scaffolds_final.agp"),
        Path(f"{prefix}.agp"),
    ]
    for c in candidates:
        if c.exists():
            return c
    globbed = sorted(glob.glob(f"{prefix}*scaffolds*final*.agp"))
    if globbed:
        return Path(globbed[0])
    raise PipelineError(
        f"could not locate yahs scaffold AGP output near prefix {prefix} "
        f"(tried {', '.join(str(c) for c in candidates)})"
    )


def step_collect_results(ctx: Ctx, run_yahs_fn) -> Step:
    def _outputs():
        outs = [ctx.results_dir / "haps.scfs.agp", ctx.results_dir / "haps.group.txt"]
        if run_yahs_fn():
            outs += [
                ctx.results_dir / "allhaps.scf.agp",
                ctx.results_dir / "allhaps.scf.hic.png",
                ctx.results_dir / "allhaps.scf.hic.pdf",
                ctx.results_dir / "allhaps.scf.fa",
            ]
        return outs

    def _run():
        ctx.results_dir.mkdir(parents=True, exist_ok=True)
        hap_prefix = ctx.hapscotch_dir / "haps"
        shutil.copy(f"{hap_prefix}.scfs.agp", ctx.results_dir / "haps.scfs.agp")
        shutil.copy(f"{hap_prefix}.group.txt", ctx.results_dir / "haps.group.txt")
        if run_yahs_fn():
            final_agp = ctx.yahs_dir / "allhaps.scf.agp"
            shutil.copy(final_agp, ctx.results_dir / "allhaps.scf.agp")
            shutil.copy(ctx.yahs_dir / "allhaps.scf.hic.png", ctx.results_dir / "allhaps.scf.hic.png")
            shutil.copy(ctx.yahs_dir / "allhaps.scf.hic.pdf", ctx.results_dir / "allhaps.scf.hic.pdf")
            final_fa = ctx.results_dir / "allhaps.scf.fa"
            run_cmd(
                [ctx.seqtools, "seq", "-o", str(final_fa),
                 str(ctx.results_dir / "allhaps.scf.agp"), str(ctx.seqfile)],
                ctx.logdir / "results.seqtools_seq.log",
            )

    return Step("collect_results", _outputs, _run)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("seqfile", help="genome sequence file (fasta[.gz]), positional")

    p.add_argument("-o", "--outdir", default=None,
                    help="output directory [HapScotch_OUT_<YYYYMMDD>]")
    p.add_argument("-t", "--threads", type=int, default=8, help="threads [8]")
    p.add_argument("-v", "--verbose", type=int, default=0, help="verbose level [0]")

    aln = p.add_argument_group("alignment inputs")
    aln.add_argument("--seq-aln", metavar="PAF",
                      help="existing self-alignment PAF; skips self_align step")
    aln.add_argument("--hic-aln", action="append", metavar="FILE", default=[],
                      help="pre-aligned HiC file (BED/PA5/BAM/BIN); repeatable")
    aln.add_argument("--hic-file", action="append", metavar="FILE", default=[],
                      help="raw HiC fastq/a file to align via hicaln.py; repeatable")

    hs = p.add_argument_group("hapscotch options")
    hs.add_argument("-g", "--busco", metavar="TSV", help="BUSCO full table")
    hs.add_argument("-p", "--ploidy", type=int, default=None, help="genome ploidy (0=auto)")
    hs.add_argument("-B", "--min-ext", type=int, default=None)
    hs.add_argument("-q", "--min-qual", type=int, default=None)
    hs.add_argument("--max-ploidy", type=int, default=None)
    hs.add_argument("--hapscotch-opt", action="append", default=[],
                     help="extra raw option passed through to hapscotch; repeatable ('hapscotch -h' for all options)")
    
    hc = p.add_argument_group("hapcure options")
    hc.add_argument("--no-contig-ec", dest="contig_ec", action="store_false", default=True,
                     help="skip hapcure contig error-correction step")
    hc.add_argument("--hapcure-opt", action="append", default=[],
                     help="extra raw option passed through to hapcure; repeatable "
                          "('hapcure -h' for all options)")
    
    yh = p.add_argument_group("yahs options")
    yh_grp = yh.add_mutually_exclusive_group()
    yh_grp.add_argument("--yahs", dest="run_yahs", action="store_true", default=None,
                         help="force-enable the YaHS branch (errors if no HiC data)")
    yh_grp.add_argument("--no-yahs", dest="run_yahs", action="store_false",
                         help="force-disable the YaHS branch even if HiC data is available")
    yh.add_argument("--yahs-opt", action="append", default=[],
                     help="extra raw option passed through to yahs; repeatable "
                          "('yahs -h' for all options)")
    
    hm = p.add_argument_group("hicmap options")
    hm.add_argument("--hictools-prepare-opt", action="append", default=[],
                     help="extra raw option passed through to hictools prepare; repeatable "
                          "('hictools prepare -h' for all options)")
    hm.add_argument("--hicmap-opt", action="append", default=[],
                     help="extra raw option passed through to hicmap.py; repeatable "
                          "('hicmap.py -h' for all options)")

    ctl = p.add_argument_group("pipeline control")
    ctl.add_argument("--resume", action="store_true",
                      help="resume a previous run in OUTDIR: skips already-completed "
                           "steps and continues from the first incomplete one, "
                           "without rerunning anything that already succeeded")
    ctl.add_argument("--force", action="store_true",
                      help="ignore existing outputs and rerun every step")
    ctl.add_argument("--force-from", metavar="STEP",
                      help="rerun this step and every step after it, by name")

    bins = p.add_argument_group(
        "tool locations",
        "bare names resolve via PATH, which automatically includes the repo root "
        "(parent of scripts/) so the compiled binaries are found without setting these",
    )
    bins.add_argument("--seqtools-bin", default="seqtools")
    bins.add_argument("--hictools-bin", default="hictools")
    bins.add_argument("--hapscotch-bin", default="hapscotch")
    bins.add_argument("--hapcure-bin", default="hapcure")
    bins.add_argument("--yahs-bin", default="yahs")
    bins.add_argument("--fastga-bin", default="FastGA")

    return p


STEP_NAMES = [
    "seq_index", "self_align", "hic_align", "hic_convert",
    "contig_ec", "hapscotch", "yahs_scaffold", "collect_results",
]

# argparse fields that don't affect what the pipeline does/produces, so they're
# excluded from the OUTDIR/data/CMD record used to validate a --force*/--resume run
_NON_CRITICAL_PARAM_KEYS = {
    "resume", "force", "force_from", "verbose", "outdir",
    "seqtools_bin", "hictools_bin", "hapscotch_bin", "hapcure_bin", "yahs_bin", "fastga_bin",
}


def critical_params(args: argparse.Namespace) -> dict:
    return {k: v for k, v in sorted(vars(args).items()) if k not in _NON_CRITICAL_PARAM_KEYS}


def diff_params(recorded: dict, current: dict) -> list:
    diffs = []
    for k in sorted(set(recorded) | set(current)):
        if recorded.get(k) != current.get(k):
            diffs.append(f"{k}: recorded={recorded.get(k)!r} vs now={current.get(k)!r}")
    return diffs


def main(argv=None) -> int:
    _add_repo_root_to_path()
    args = build_arg_parser().parse_args(argv)

    outdir = Path(args.outdir) if args.outdir else Path(
        f"HapScotch_OUT_{__import__('datetime').date.today():%Y%m%d}"
    )

    resuming = bool(args.resume or args.force or args.force_from)
    outdir_existed = outdir.exists()

    if outdir_existed and not resuming:
        _log(f"error: output directory '{outdir}' already exists. Use --resume to "
             f"continue from where it stopped, --force to rerun everything, or "
             f"--force-from STEP to rerun from a specific step.")
        return 1

    if resuming and not outdir_existed:
        _log(f"error: --resume/--force/--force-from given but output directory "
             f"'{outdir}' does not exist - nothing to resume.")
        return 1

    if args.force_from and args.force_from not in STEP_NAMES:
        _log(f"error: unknown --force-from step '{args.force_from}', "
             f"expected one of: {', '.join(STEP_NAMES)}")
        return 1

    ctx = Ctx(
        seqfile=Path(args.seqfile), outdir=outdir, threads=args.threads, verbose=args.verbose,
        seqtools=args.seqtools_bin, hictools=args.hictools_bin,
        hapscotch_bin=args.hapscotch_bin, hapcure_bin=args.hapcure_bin,
        yahs_bin=args.yahs_bin,
    )
    ctx.ensure_dirs()

    cmd_file = ctx.datadir / "CMD"
    current_params = critical_params(args)
    if resuming:
        try:
            recorded_params = json.loads(cmd_file.read_text())
        except (OSError, json.JSONDecodeError):
            _log(f"error: cannot resume - no valid parameter record found at {cmd_file} "
                 f"(remove {outdir} to start a fresh run instead)")
            return 1
        mismatches = diff_params(recorded_params, current_params)
        if mismatches:
            _log("error: cannot resume - parameters differ from the run that created "
                 f"{outdir}:")
            for line in mismatches:
                _log(f"  {line}")
            return 1
    else:
        cmd_file.write_text(json.dumps(current_params, indent=2, sort_keys=True) + "\n")

    # --- wire up steps (later steps consume earlier steps' resolved-path callbacks) ---
    step_idx = step_seq_index(ctx)
    step_sa, seqaln_fn = step_self_align(ctx, args)
    step_ha, resolved_hicalns_fn = step_hic_align(ctx, args)
    step_hc, resolved_hicbin_fn, will_have_hic_fn = step_hic_convert(ctx, resolved_hicalns_fn)
    step_ec, resolved_agpec_fn = step_contig_ec(ctx, args, resolved_hicbin_fn, will_have_hic_fn)

    def run_yahs_fn() -> bool:
        if args.run_yahs is False:
            return False
        have_hic = will_have_hic_fn()
        if args.run_yahs is True and not have_hic:
            raise PipelineError("--yahs was forced but no HiC data is available "
                                 "(need --hic-aln/--hic-file)")
        return bool(args.run_yahs) if args.run_yahs is not None else have_hic

    step_hs = step_hapscotch(ctx, args, seqaln_fn, resolved_hicbin_fn, resolved_agpec_fn,
                              run_yahs_fn)
    step_yh = step_yahs_scaffold(ctx, args, resolved_hicbin_fn, run_yahs_fn)
    step_res = step_collect_results(ctx, run_yahs_fn)

    steps = [step_idx, step_sa, step_ha, step_hc, step_ec, step_hs, step_yh, step_res]

    force_from_seen = False
    try:
        preflight_check(ctx, steps)
        for step in steps:
            if args.force_from and step.name == args.force_from:
                force_from_seen = True
            run_step(step, force=args.force, force_from=force_from_seen)
    except PipelineError as e:
        _log(f"ERROR: {e}")
        return 1

    _log(f"pipeline complete. Results in {ctx.results_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
