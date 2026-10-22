#!/usr/bin/env python3
"""
run_pipeline.py - end-to-end driver for the HapScotch haplotyping pipeline.

Stages (each is independently resumable - see Checkpointing below):
  1. seq_index       seqtools idx                   -> genome sequence index (.idx)
  2. self_align      selfaln.py                     -> self-alignment (.1aln)
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
already holds a previous run (OUTDIR/data/CMD exists), one of --resume/--force/
--force-from is required (otherwise the pipeline refuses to run); conversely, those
options require such a previous run. An existing OUTDIR without a previous run
is used as is for a new run. The parameters of the run that created OUTDIR are recorded
in OUTDIR/data/CMD; --resume/--force/--force-from are refused if the current
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
from datetime import datetime
from pathlib import Path
from typing import Callable, Optional

# --------------------------------------------------------------------------- #
# errors / logging
# --------------------------------------------------------------------------- #

class PipelineError(RuntimeError):
    """Raised when a step fails or its expected outputs are missing after running."""


def _log(msg: str) -> None:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[run_pipeline] [{ts}] {msg}", file=sys.stderr, flush=True)


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


def _resolve_fastga_tool(tool: str, fastga_bin: str) -> str:
    """Resolve a FastGA companion tool (ALNchain, ALNplot): a bare name is looked up
    next to the FastGA executable first, then on PATH."""
    if os.sep not in tool:
        sibling = Path(fastga_bin).parent / tool
        if Path(fastga_bin).parent != Path(".") and os.access(sibling, os.X_OK):
            return str(sibling.resolve())
    return _resolve_exe(tool)


def _resolve_exe(tool: str) -> str:
    """Resolve tool name or path to its real absolute path if it exists, else keep unchanged."""
    found = shutil.which(tool)
    if found:
        return str(Path(found).resolve())
    p = Path(tool)
    if p.exists():
        return str(p.resolve())
    return tool


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
    fastga_bin: str = "FastGA"
    bwamem2_bin: str = "bwa-mem2"
    minibwa_bin: str = "minibwa"
    hic_aligner: str = "bwamem2"
    seq_aligner: str = "FastGA"
    samtools_bin: str = "samtools"
    alnsplit_bin: str = "alnsplit"
    alnchain_bin: str = "ALNchain"
    alnplot_bin: str = "ALNplot"

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
    def hicaln_dir(self) -> Path:
        return self.datadir / "hicalns"

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


def link_ref_genome(ctx: Ctx) -> Path:
    """Symlink the input genome into datadir/ref.fa[.gz] and return that path,
    for use as the reference from here on instead of ctx.seqfile directly."""
    target = ctx.seqfile.resolve()
    ref_link = ctx.datadir / ("ref.fa.gz" if target.name.endswith(".gz") else "ref.fa")
    if ref_link.is_symlink():
        if ref_link.resolve() != target:
            raise PipelineError(
                f"{ref_link} already exists and points elsewhere ({ref_link.resolve()}), "
                f"expected it to point to {target} - remove it or use a fresh OUTDIR"
            )
    elif ref_link.exists():
        raise PipelineError(f"{ref_link} already exists and is not a symlink - "
                             f"remove it or use a fresh OUTDIR")
    else:
        ref_link.symlink_to(target)
    return ref_link


# --------------------------------------------------------------------------- #
# command execution
# --------------------------------------------------------------------------- #

def run_cmd(cmd, log_path: Path, *, shell: bool = False, cwd: Optional[Path] = None,
            env: Optional[dict] = None) -> None:
    """Run `cmd`, tee-ing stdout/stderr to log_path. `env` adds to (does not replace)
    the current environment for this command only."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    printable = cmd if isinstance(cmd, str) else " ".join(shlex.quote(str(c)) for c in cmd)
    _log(f"running: {printable}")
    _log(f"  log: {log_path}")
    with open(log_path, "w") as lf:
        lf.write(f"# CMD: {printable}\n")
        lf.flush()
        try:
            proc = subprocess.run(
                cmd, shell=shell, cwd=cwd, stdout=lf, stderr=subprocess.STDOUT,
                env=None if env is None else {**os.environ, **env},
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


def _parse_opts(opt_list: list) -> list:
    """Splits space-delimited option strings into individual command-line arguments."""
    res = []
    for opt_str in opt_list:
        res.extend(shlex.split(opt_str))
    return res


def require_outputs(outputs) -> None:
    missing = [str(o) for o in outputs if not Path(o).exists()]
    if missing:
        raise PipelineError(
            "step completed but expected output(s) missing: " + ", ".join(missing)
        )


# --------------------------------------------------------------------------- #
# preflight - check tool availability before running
# --------------------------------------------------------------------------- #

def _check_binary(name: str, target: str) -> Optional[str]:
    if shutil.which(target) is None:
        return f"{name}: '{target}' not found on PATH (or not executable)"
    return None


# ALNplot needs one of these to convert its EPS output to PDF
PDF_TOOLS = ("pstopdf", "epstopdf", "ps2pdf", "eps2pdf")


def _check_pdf_tool(name: str, target) -> Optional[str]:
    if not any(shutil.which(t) for t in PDF_TOOLS):
        return f"{name}: none of {'/'.join(PDF_TOOLS)} found on PATH (needed by ALNplot)"
    return None


def _check_script(name: str, target: Path) -> Optional[str]:
    if not target.is_file():
        return f"{name}: expected script not found: {target}"
    return None


def preflight_check(ctx: Ctx, steps: list, aln_plot: bool = False) -> None:
    """Verify every external tool/script needed by applicable steps can be found."""
    aligner_bin = ctx.bwamem2_bin if ctx.hic_aligner == "bwamem2" else ctx.minibwa_bin

    checks_by_step = {
        "seq_index":       [("seqtools", ctx.seqtools)],
        "self_align":      [("selfaln.py", SELFALN_PY), (ctx.fastga_bin, ctx.fastga_bin)],
        "hic_align":       [("hicaln.py", HICALN_PY),
                            (ctx.hic_aligner, aligner_bin),
                            ("samtools", ctx.samtools_bin)],
        "hic_convert":     [("hictools", ctx.hictools)],
        "contig_ec":       [("hapcure", ctx.hapcure_bin)],
        "hapscotch":       [("hapscotch", ctx.hapscotch_bin)],
        "yahs_scaffold":   [("seqtools", ctx.seqtools), ("yahs", ctx.yahs_bin),
                            ("hictools", ctx.hictools), ("hicmap.py", HICMAP_PY)],
        "collect_results": [("seqtools", ctx.seqtools)],
    }
    if aln_plot:
        checks_by_step["collect_results"] += [
            ("alnsplit", ctx.alnsplit_bin), ("ALNchain", ctx.alnchain_bin),
            ("ALNplot", ctx.alnplot_bin), ("[e]ps[to|2]pdf", PDF_TOOLS)]

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
            if isinstance(target, Path):
                err = _check_script(name, target)
            elif isinstance(target, tuple):
                err = _check_pdf_tool(name, target)
            else:
                err = _check_binary(name, target)
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
    outputs: Callable[[], list]
    run: Callable[[], None]
    applicable: Callable[[], bool] = field(default=lambda: True)


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


def step_self_align(ctx: Ctx, args) -> tuple:
    out = ctx.datadir / "seqaln.1aln"

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
    prefixes = [ctx.hicaln_dir / f"hicaln.{i}" for i in range(len(hicfiles))]
    generated = [Path(f"{prefix}.bam") for prefix in prefixes]

    def _applicable() -> bool:
        return len(hicfiles) > 0

    def _run():
        ctx.hicaln_dir.mkdir(parents=True, exist_ok=True)
        for src, prefix in zip(hicfiles, prefixes):
            run_cmd(
                [sys.executable, str(HICALN_PY), str(ctx.seqfile), str(src),
                 "-o", str(prefix), "-t", str(ctx.threads),
                 "-a", args.hic_aligner,
                 "--bwamem2-bin", args.bwamem2_bin,
                 "--minibwa-bin", args.minibwa_bin,
                 "--samtools-bin", args.samtools_bin],
                ctx.logdir / f"hic_align.{prefix.name}.log",
            )

    def _resolved_hicalns() -> list:
        return list(args.hic_aln or []) + generated

    return Step("hic_align", lambda: list(generated), _run, applicable=_applicable), _resolved_hicalns


def step_hic_convert(ctx: Ctx, resolved_hicalns_fn) -> tuple:
    out = ctx.datadir / "hicaln.bin"

    def _is_passthrough() -> bool:
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
        if _is_passthrough():
            return Path(hicalns[0])
        return None

    def _will_have_hic() -> bool:
        return len(resolved_hicalns_fn()) > 0

    return Step("hic_convert", _outputs, _run,
                applicable=_applicable), _resolved_hicbin, _will_have_hic


def step_contig_ec(ctx: Ctx, args, resolved_hicbin_fn, will_have_hic_fn, seqaln_fn) -> tuple:
    out = ctx.hapcure_dir / "ctg.ec.agp"

    def _applicable() -> bool:
        return args.agp is None and args.contig_ec and will_have_hic_fn()

    def _run():
        ctx.hapcure_dir.mkdir(parents=True, exist_ok=True)
        run_cmd(
            [ctx.hapcure_bin, "--hic-tmpdir", str(ctx.hapcure_dir)] + _parse_opts(args.hapcure_opt) +
            ["-a", str(seqaln_fn())] +
            ["-t", str(ctx.threads)] +
            ["-o", str(ctx.hapcure_dir / "ctg")] +
            [str(ctx.idxfile), str(resolved_hicbin_fn())],
            ctx.logdir / "contig_ec.log",
        )

    def _resolved_agpec() -> Optional[Path]:
        if args.agp is not None:
            return Path(args.agp)
        return out if _applicable() else None

    return Step("contig_ec", lambda: [out], _run, applicable=_applicable), _resolved_agpec


def step_hapscotch(ctx: Ctx, args, seqaln_fn, resolved_hicbin_fn, resolved_agpec_fn,
                    run_yahs_fn) -> Step:
    prefix = ctx.hapscotch_dir / "haps"
    outs = [Path(f"{prefix}.grp.agp"),
            Path(f"{prefix}.grp.txt"),
            Path(f"{prefix}.ploidy")]

    def _run():
        ctx.hapscotch_dir.mkdir(parents=True, exist_ok=True)
        opts = ["--hic-tmpdir", str(ctx.hapscotch_dir)]
        opts += _parse_opts(args.hapscotch_opt)
        if args.ploidy is not None:
            opts += ["-p", str(args.ploidy)]
        opts += ["-t", str(ctx.threads)]
        agpec = resolved_agpec_fn()
        if agpec is not None:
            opts += ["-a", str(agpec)]
        hicbin = resolved_hicbin_fn()
        if hicbin is not None:
            opts += ["-c", str(hicbin), "--file-type", "BIN"]
        if run_yahs_fn():
            opts += ["-Y"]
        run_cmd(
            [ctx.hapscotch_bin] + opts + ["-o", str(prefix), str(ctx.idxfile), str(seqaln_fn())],
            ctx.logdir / "hapscotch.log",
        )
        
    return Step("hapscotch", lambda: outs, _run)


def step_yahs_scaffold(ctx: Ctx, args, resolved_hicbin_fn, run_yahs_fn) -> Step:
    hap_prefix = ctx.hapscotch_dir / "haps"
    yahs_prefix = ctx.yahs_dir / "bbseq"
    final_agp = f"{yahs_prefix}_scaffolds_final.agp"

    def _applicable() -> bool:
        return run_yahs_fn()

    def _run():
        ctx.yahs_dir.mkdir(parents=True, exist_ok=True)
        bbseq_agp = Path(f"{hap_prefix}.bbseq.agp")
        bbscf_agp = Path(f"{hap_prefix}.bbscf.agp")
        bbscf_hic_bin = Path(f"{hap_prefix}.bbseq-hic.bin")
        seq_hic_bin = resolved_hicbin_fn()
        for needed in (bbseq_agp, bbscf_agp, bbscf_hic_bin, seq_hic_bin):
            if not needed.exists():
                raise PipelineError(
                    f"YaHS branch requested but required hapscotch output missing: {needed}"
                )

        bbseq_fa = ctx.yahs_dir / "haps.bbseq.fa.gz"
        run_cmd(
            [ctx.seqtools, "seq", "-t", str(args.threads), "-a", "-o", str(bbseq_fa), 
             str(ctx.seqfile), str(bbseq_agp)],
             ctx.logdir / "yahs.seqtools_bbseq.log",
        )

        bbseq_idx = f"{bbseq_fa}.fai"
        run_cmd(
            [ctx.seqtools, "idx", "-o", str(bbseq_idx), str(bbseq_fa)],
            ctx.logdir / "yahs.seqtools_bbidx.log",
        )

        run_cmd(
            [ctx.yahs_bin] + _parse_opts(args.yahs_opt) + 
            ["-a", str(bbscf_agp), "--no-contig-ec", "--no-scaffold-ec",
             "-o", str(yahs_prefix), str(bbseq_fa), str(bbscf_hic_bin)],
            ctx.logdir / "yahs.log",
        )

    return Step("yahs_scaffold", lambda: [final_agp], _run, applicable=_applicable)


def _find_yahs_scaffolds_agp(prefix: Path) -> Path:
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
        f"could not locate yahs scaffold AGP output near prefix {prefix}"
    )


def _write_seq_order(scf_agp: Path, out: Path) -> int:
    """ALNplot sequence order from the W lines of a scaffold AGP: '@<name>_<beg>_<end><ori>',
    matching the piece names written by 'alnsplit -r'."""
    n = 0
    with open(scf_agp) as f, open(out, "w") as fo:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[4] != "W":
                continue
            ori = "-" if c[8] == "-" else "+"
            fo.write(f"@{c[5]}_{c[6]}_{c[7]}{ori}\n")
            n += 1
    return n


def _write_identity_agp(idxfile: Path, out: Path) -> None:
    """A break AGP with each sequence as one whole piece, for when there is no
    error-correction AGP (alnsplit then only renames sequences to name_1_len)."""
    with open(idxfile) as f, open(out, "w") as fo:
        for line in f:
            c = line.split()
            if len(c) < 2:
                continue
            fo.write(f"{c[0]}\t1\t{c[1]}\t1\tW\t{c[0]}\t1\t{c[1]}\t+\n")


def make_aln_plot(ctx: Ctx, args, seqaln: Path, agpec: Optional[Path], scf_agp: Path, pdf: Path) -> None:
    """Plot the self-alignment against the final scaffolds:
    1. alnsplit -r   map the alignments onto the error-corrected pieces, named name_beg_end
    2. ALNchain      keep the best chains, otherwise there are too many alignments to plot
    3. ALNplot       plot in scaffold order, given by the W lines of the scaffold AGP
    the sequence order and the split alignment file are removed afterwards"""
    seq_order = ctx.datadir / "seq.order"
    split_aln = ctx.datadir / "seqaln_ec.1aln"
    chain_aln = ctx.datadir / "seqaln_ec.chain.1aln"
    id_agp = ctx.datadir / "seq.id.agp"

    if _write_seq_order(scf_agp, seq_order) == 0:
        raise PipelineError(f"no sequence components found in scaffold AGP: {scf_agp}")
    if agpec is None:
        _write_identity_agp(ctx.idxfile, id_agp)
    try:
        run_cmd(
            [ctx.alnsplit_bin, "-r", "-o", str(split_aln),
             str(ctx.idxfile), str(agpec if agpec is not None else id_agp), str(seqaln)],
            ctx.logdir / "res.alnsplit.log",
        )
        # FastGA tools take option values attached, e.g. -o<file>
        run_cmd(
            [ctx.alnchain_bin, f"-o{chain_aln}", str(split_aln)],
            ctx.logdir / "res.alnchain.log",
        )
        # ALNplot writes <pdf>.eps then converts it with [e]ps[to|2]pdf, and ps2pdf puts
        # the PDF in the working directory - so run it from the PDF's own directory.
        # ps2pdf also ignores the EPS bounding box and uses a default (Letter/A4) page,
        # clipping the square plot; -dEPSCrop (via GS_OPTIONS, read by Ghostscript)
        # crops the page to the bounding box as epstopdf does
        gs_opts = " ".join(o for o in (os.environ.get("GS_OPTIONS", ""), "-dEPSCrop") if o)
        run_cmd(
            [ctx.alnplot_bin, "-L", "-G", "-n0", f"-T{ctx.threads}", f"-p:{pdf.name}",
             str(chain_aln.resolve()), str(seq_order.resolve()), str(seq_order.resolve())],
            ctx.logdir / "res.alnplot.log", cwd=pdf.parent,
            env={"GS_OPTIONS": gs_opts},
        )
        pdf.with_suffix(".eps").unlink(missing_ok=True)
    finally:
        for f in (seq_order, split_aln, id_agp):
            f.unlink(missing_ok=True)


def step_collect_results(ctx: Ctx, args, resolved_agpec_fn, resolved_hicbin_fn, run_yahs_fn,
                          seqaln_fn) -> Step:
    def _aln_plot() -> bool:
        return args.aln_plot and run_yahs_fn()

    def _outputs():
        outs = [ctx.results_dir / "haps.cnt.txt",
                ctx.results_dir / "haps.grp.agp", 
                ctx.results_dir / "haps.grp.txt"]
        if resolved_agpec_fn():
            outs += [ctx.results_dir / "haps.ctg-ec.agp"]
        if run_yahs_fn():
            outs += [
                ctx.results_dir / "haps-all.scf.agp",
                ctx.results_dir / "haps-all.scf.hic.png",
                ctx.results_dir / "haps-all.scf.hic.pdf",
            ]
        if _aln_plot():
            outs += [ctx.results_dir / "haps-all.scf.aln.pdf"]
        return outs

    def _run():
        ctx.results_dir.mkdir(parents=True, exist_ok=True)
        hap_prefix = ctx.hapscotch_dir / "haps"
        shutil.copy(f"{hap_prefix}.grp.agp", ctx.results_dir / "haps.grp.agp")
        shutil.copy(f"{hap_prefix}.grp.txt", ctx.results_dir / "haps.grp.txt")
        shutil.copy(f"{hap_prefix}.ploidy",  ctx.results_dir / "haps.cnt.txt")

        ploidy = 0
        with open(ctx.results_dir / "haps.cnt.txt") as f:
            ploidy = max((int(line.split()[1]) for line in f), default=0)

        grp_file = str(ctx.results_dir / "haps.grp.txt")
        for hap in range(1, ploidy + 1):
            hap_agp = ctx.results_dir / f"haps-{hap}.ctg.agp"
            hap_fa  = ctx.results_dir / f"haps-{hap}.ctg.fa.gz"
            with open(hap_agp, "w") as out:
                with open(grp_file, "r") as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("#"):
                            continue
                        parts = line.split()
                        if len(parts) < 6 or int(parts[4]) != hap:
                            continue
                        s_id  = parts[0]
                        s_beg = int(parts[1])
                        s_end = int(parts[2])
                        comp_len = s_end - s_beg
                        agp_line = f"{s_id}_{s_beg+1}_{s_end}\t1\t{comp_len}\t1\tW\t{s_id}\t{s_beg+1}\t{s_end}\t+\n"
                        out.write(agp_line)
            run_cmd(
                [ctx.seqtools, "seq", "-t", str(args.threads), "-a", "-o", str(hap_fa), 
                 str(ctx.seqfile), str(hap_agp)],
                 ctx.logdir / f"res.seqtools_ctg_hap{hap}_fa.log",
            )

        agpec = resolved_agpec_fn()
        if agpec is not None:
            shutil.copy(str(agpec),  ctx.results_dir / "haps.ctg-ec.agp")    
        
        if run_yahs_fn():
            yahs_agp = _find_yahs_scaffolds_agp(ctx.yahs_dir / "bbseq")
            bbpos_txt = ctx.hapscotch_dir / "haps.bbpos.txt"
            
            scfall_agp = ctx.results_dir / "haps-all.scf.agp"
            run_cmd(
                [ctx.seqtools, "hap", "-o", str(scfall_agp), 
                 str(yahs_agp), str(bbpos_txt)],
                ctx.logdir / "res.seqtools_scf_all_agp.log",
            )

            for hap in range(1, ploidy + 1):
                hap_agp = ctx.results_dir / f"haps-{hap}.scf.agp"
                hap_fa  = ctx.results_dir / f"haps-{hap}.scf.fa.gz"
                run_cmd(
                    [ctx.seqtools, "hap", "-p", str(hap), "-o", str(hap_agp), 
                     str(yahs_agp), str(bbpos_txt)],
                    ctx.logdir / f"res.seqtools_scf_hap{hap}_agp.log",
                )
                run_cmd(
                    [ctx.seqtools, "seq", "-t", str(args.threads), "-a", "-o", str(hap_fa),
                     str(ctx.seqfile), str(hap_agp)],
                    ctx.logdir / f"res.seqtools_scf_hap{hap}_fa.log",
                )

            hic_txt = ctx.results_dir / "haps-all.scf.hic.txt"
            run_cmd(
                [ctx.hictools, "prepare", "--hic-tmpdir", str(ctx.results_dir), 
                 "-a", str(scfall_agp), "-n", "3000", "-o", str(hic_txt)] + 
                [str(resolved_hicbin_fn()), str(ctx.idxfile)],
                ctx.logdir / "res.hictools_prepare.log",
            )

            hic_png = ctx.results_dir / "haps-all.scf.hic.png"
            hic_pdf = ctx.results_dir / "haps-all.scf.hic.pdf"
            run_cmd(
                [sys.executable, str(HICMAP_PY), 
                 "--png", str(hic_png), "--pdf", str(hic_pdf), 
                 str(hic_txt)],
                ctx.logdir / "res.hicmap.log",
            )

            if _aln_plot():
                make_aln_plot(ctx, args, Path(seqaln_fn()), agpec, scfall_agp,
                              ctx.results_dir / "haps-all.scf.aln.pdf")

    return Step("collect_results", _outputs, _run)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("seqfile", help="genome sequence file (fasta[.gz]), positional")

    p.add_argument("-o", "--outdir", default=None, help="output directory [HapScotch_OUT_<YYYYMMDD>]")
    p.add_argument("-a", "--agp", metavar="AGP", 
                   help="existing error-correction AGP; skips error correction")
    p.add_argument("-p", "--ploidy", type=int, default=0, 
                   help="ploidy number; skips ploidy estimation")
    p.add_argument("-t", "--threads", type=int, default=8, help="threads [8]")
    p.add_argument("-v", "--verbose", type=int, default=0, help="verbose level [0]")

    aln = p.add_argument_group("alignment inputs")
    aln.add_argument("--seq-aln", metavar="FILE",
                      help="existing self-alignment (.1aln recommended, or PAF[.gz]); skips self_align step")
    aln.add_argument("--hic-aln", action="append", metavar="FILE", default=[],
                      help="pre-aligned HiC file (BED/PA5/BAM/BIN); repeatable")
    aln.add_argument("--hic-file", action="append", metavar="FILE", default=[],
                      help="raw HiC fastq/a file to align via hicaln.py; repeatable")
    
    ap = p.add_argument_group("aligner program")
    ap.add_argument("--seq-aligner", choices=["FastGA"], default="FastGA",
                      help="aligner program for sequence comparison [FastGA]")
    ap.add_argument("--hic-aligner", choices=["bwamem2", "minibwa"], default="bwamem2",
                      help="aligner program for raw HiC files [bwamem2]")

    hs = p.add_argument_group("hapscotch options")
    hs.add_argument("--hapscotch-opt", action="append", default=[],
                     help="options passed through to hapscotch; repeatable ('hapscotch -h' for all options)")
    
    hc = p.add_argument_group("hapcure options")
    hc.add_argument("--no-contig-ec", dest="contig_ec", action="store_false", default=True,
                     help="skip hapcure contig error-correction step")
    hc.add_argument("--hapcure-opt", action="append", default=[],
                     help="options passed through to hapcure; repeatable ('hapcure -h' for all options)")
    
    yh = p.add_argument_group("yahs options")
    yh.add_argument("--no-yahs", dest="run_yahs", action="store_false", default=True,
                         help="force-disable the YaHS branch even if HiC data is available")
    yh.add_argument("--yahs-opt", action="append", default=[],
                     help="options passed through to yahs; repeatable ('yahs -h' for all options)")

    rs = p.add_argument_group("result options")
    rs.add_argument("--no-aln-plot", dest="aln_plot", action="store_false", default=True,
                     help="skip the sequence alignment plot of the scaffolds (needs alnsplit, "
                          "and ALNchain/ALNplot from FastGA)")

    ctl = p.add_argument_group("pipeline control")
    ctl.add_argument("--resume", action="store_true",
                      help="resume a previous run in OUTDIR: skips already-completed steps")
    ctl.add_argument("--force", action="store_true",
                      help="ignore existing outputs and rerun every step")
    ctl.add_argument("--force-from", metavar="STEP",
                      help="rerun this step and every step after it, by name")

    bins = p.add_argument_group(
        "tool locations",
        "bare names resolve via PATH (parent of scripts/ automatically included)",
    )
    bins.add_argument("--seqtools-bin", default="seqtools")
    bins.add_argument("--hictools-bin", default="hictools")
    bins.add_argument("--hapscotch-bin", default="hapscotch")
    bins.add_argument("--hapcure-bin", default="hapcure")
    bins.add_argument("--yahs-bin", default="yahs")
    bins.add_argument("--fastga-bin", default="FastGA")
    bins.add_argument("--bwamem2-bin", default="bwa-mem2")
    bins.add_argument("--minibwa-bin", default="minibwa")
    bins.add_argument("--samtools-bin", default="samtools")
    bins.add_argument("--alnsplit-bin", default="alnsplit")
    bins.add_argument("--alnchain-bin", default="ALNchain")
    bins.add_argument("--alnplot-bin", default="ALNplot")

    return p


STEP_NAMES = [
    "seq_index", "self_align", "hic_align", "hic_convert",
    "contig_ec", "hapscotch", "yahs_scaffold", "collect_results",
]

_NON_CRITICAL_PARAM_KEYS = {
    "resume", "force", "force_from", "verbose", "outdir",
    "seqtools_bin", "hictools_bin", "hapscotch_bin", "hapcure_bin", "yahs_bin", "fastga_bin",
    "bwamem2_bin", "minibwa_bin", "samtools_bin",
    "alnsplit_bin", "alnchain_bin", "alnplot_bin",
    # result-only option, added later; not recorded by older runs
    "aln_plot",
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

    for k in ("seqtools_bin", "hictools_bin", "hapscotch_bin", "hapcure_bin",
              "yahs_bin", "fastga_bin", "bwamem2_bin", "minibwa_bin", "samtools_bin",
              "alnsplit_bin"):
        setattr(args, k, _resolve_exe(getattr(args, k)))
    # FastGA companion tools usually live next to the FastGA executable
    for k in ("alnchain_bin", "alnplot_bin"):
        setattr(args, k, _resolve_fastga_tool(getattr(args, k), args.fastga_bin))

    outdir = Path(args.outdir) if args.outdir else Path(
        f"HapScotch_OUT_{datetime.now():%Y%m%d}"
    )

    resuming = bool(args.resume or args.force or args.force_from)

    # an existing OUTDIR without a parameter record (e.g. pre-created by a job
    # runner) has never been used by a pipeline run, so it is treated as a new one
    outdir_existed = (outdir / "data" / "CMD").exists()
    if outdir.exists() and not outdir_existed and not resuming:
        _log(f"output directory '{outdir}' exists but has no previous run record, "
             f"starting a new run in it")

    if outdir_existed and not resuming:
        _log(f"error: output directory '{outdir}' already exists. Use --resume to "
             f"continue from where it stopped, --force to rerun everything, or "
             f"--force-from STEP to rerun a specific step.")
        return 1

    if resuming and not outdir_existed:
        _log(f"error: --resume/--force/--force-from given but output directory "
             f"'{outdir}' does not exist - nothing to resume.")
        return 1

    if args.force_from and args.force_from not in STEP_NAMES:
        _log(f"error: unknown --force-from step '{args.force_from}', "
             f"expected one of: {', '.join(STEP_NAMES)}")
        return 1

    if args.agp is not None and not Path(args.agp).is_file():
        _log(f"error: error-correction AGP not found or not a regular file: {args.agp}")
        return 1

    ctx = Ctx(
        seqfile=Path(args.seqfile), outdir=outdir, threads=args.threads, verbose=args.verbose,
        seqtools=args.seqtools_bin, hictools=args.hictools_bin,
        hapscotch_bin=args.hapscotch_bin, hapcure_bin=args.hapcure_bin,
        yahs_bin=args.yahs_bin, fastga_bin=args.fastga_bin,
        bwamem2_bin=args.bwamem2_bin, minibwa_bin=args.minibwa_bin,
        seq_aligner=args.seq_aligner, hic_aligner=args.hic_aligner, samtools_bin=args.samtools_bin,
        alnsplit_bin=args.alnsplit_bin, alnchain_bin=args.alnchain_bin, alnplot_bin=args.alnplot_bin,
    )
    ctx.ensure_dirs()

    try:
        ctx.seqfile = link_ref_genome(ctx)
    except PipelineError as e:
        _log(f"ERROR: {e}")
        return 1

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

    step_idx = step_seq_index(ctx)
    step_sa, seqaln_fn = step_self_align(ctx, args)
    step_ha, resolved_hicalns_fn = step_hic_align(ctx, args)
    step_hc, resolved_hicbin_fn, will_have_hic_fn = step_hic_convert(ctx, resolved_hicalns_fn)
    step_ec, resolved_agpec_fn = step_contig_ec(ctx, args, resolved_hicbin_fn, will_have_hic_fn, seqaln_fn)

    def run_yahs_fn() -> bool:
        if args.run_yahs is False:
            return False
        return will_have_hic_fn()

    step_hs = step_hapscotch(ctx, args, seqaln_fn, resolved_hicbin_fn, resolved_agpec_fn, run_yahs_fn)
    step_yh = step_yahs_scaffold(ctx, args, resolved_hicbin_fn, run_yahs_fn)
    step_res = step_collect_results(ctx, args, resolved_agpec_fn, resolved_hicbin_fn, run_yahs_fn,
                                     seqaln_fn)

    steps = [step_idx, step_sa, step_ha, step_hc, step_ec, step_hs, step_yh, step_res]

    force_from_seen = False
    try:
        preflight_check(ctx, steps, aln_plot=args.aln_plot and run_yahs_fn())
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
