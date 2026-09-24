# HapScotch - haplotype resolved polyploid genome scaffolding 
![Language](https://img.shields.io/badge/language-C-blue.svg)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey.svg)

## Overview

HapScotch phases and scaffolds haplotypes in polyploid/multi-haplotype genome assemblies from all-vs-all sequence self-alignments, optionally refined with Hi-C data (contig error-correction and Hi-C guided rescaffolding via [YaHS](https://github.com/c-zhou/yahs)).

The repository builds five C programs, plus a Python script that runs the full pipeline end to end.

| Program | Purpose |
|---|---|
| `hapscotch` | Core haplotype binning/phasing/scaffolding |
| `hapcure` | Hi-C contig error-correction |
| `hictools` | Hi-C alignment file processing |
| `seqtools` | Sequence/AGP utilities |
| `hapcount` | Standalone ploidy estimator |
| `scripts/run_pipeline.py` | Driver script that runs the tools end to end |

## Table of Contents

- [Installation](#installation)
  - [From source](#from-source)
  - [Conda](#conda)
  - [Docker / Singularity / Apptainer](#docker--singularity--apptainer)
- [Running the pipeline](#running-the-pipeline)
- [Running step by step](#running-step-by-step)
- [Other tools](#other-tools)

## Installation

You need a C compiler, GNU Make, CMake, and the zlib/pthread development headers.

### From source

```bash
git clone --recurse-submodules https://github.com/c-zhou/hapscotch.git
cd hapscotch
make
make install   # optionally installs all executables in ~/bin/
```

If you already cloned without `--recurse-submodules`, run `git submodule update --init --recursive` first — this fetches the bundled [HiGHS](https://github.com/ERGO-Code/HiGHS) ILP solver used for Hi-C phasing.

This builds `hapscotch`, `hapcure`, `hictools`, `seqtools`, and `hapcount`. It does **not** install the external pipeline tools (FastGA / minibwa / samtools / YaHS) — see the options below if you want everything set up in one step.

### Conda

**Option A — build HapScotch as a conda package.** [`recipes/conda/`](recipes/conda/) is a `conda-build` recipe that pulls in FastGA/minibwa/samtools/YaHS as run dependencies, so a single build gives you everything:

```bash
# create a conda environment with all tools
conda build recipes/conda/
conda create -n hapscotch --use-local hapscotch

# activate the environment and run the pipeline
conda activate hapscotch
run_pipeline.py --help
```

**Option B — set up dependencies only, build HapScotch yourself.** [`recipes/environment.yml`](recipes/environment.yml) sets up the build toolchain plus FastGA/minibwa/samtools/YaHS and the Python packages used by `hicmap.py`:

```bash
# create a conda environment with all dependencies
conda env create -f recipes/environment.yml
conda activate hapscotch

# build hapscotch
make
make install   # optionally installs all executables in ~/bin/

# run the pipeline
python3 scripts/run_pipeline.py --help
```

### Docker / Singularity / Apptainer

**Docker** — [`recipes/Dockerfile`](recipes/Dockerfile) builds an image with the full pipeline (all five binaries, `scripts/*.py`, and FastGA/minibwa/samtools/YaHS):

```bash
docker build -t hapscotch -f recipes/Dockerfile .
docker run --rm -v "$PWD":/data -w /data hapscotch run_pipeline.py --help
```

**Singularity/Apptainer** — [`recipes/Singularity.def`](recipes/Singularity.def) is the equivalent definition, built from the same conda-forge/bioconda base:

```bash
# or replace apptainer with singularity
apptainer build hapscotch.sif recipes/Singularity.def
apptainer exec hapscotch.sif run_pipeline.py --help
```

## Running the pipeline

The `scripts/run_pipeline.py` script runs the whole pipeline and is resumable:

| Flag | Effect |
|---|---|
| `--resume` | Resume from the failed step |
| `--force-from STEP` | Force a rerun from a given step |
| `--force` | Force a rerun from scratch |

**Prerequisites:**
- The five binaries above, either on `PATH` or pointed to with `--seqtools-bin` / `--hictools-bin` / `--hapscotch-bin` / `--hapcure-bin`
- [YaHS](https://github.com/c-zhou/yahs) (`--yahs-bin`) — only needed for the Hi-C rescaffolding branch
- [FastGA](https://github.com/thegenemyers/FASTGA) (`--fastga-bin`) — needed unless you already have self-alignments between sequences
- [minibwa](https://github.com/lh3/minibwa) (`--minibwa-bin`) and [samtools](https://github.com/samtools/samtools) (`--samtools-bin`) — needed unless you already have Hi-C alignments

**Basic usage (haploid / no Hi-C):**

```bash
python3 scripts/run_pipeline.py -o OUTDIR genome.fa
```

**With Hi-C data** (enables contig error-correction and, by default, the YaHS rescaffolding branch):

```bash
python3 scripts/run_pipeline.py -o OUTDIR --hic-aln hic.bam genome.fa
python3 scripts/run_pipeline.py -o OUTDIR --hic-file r1.fq.gz,r2.fq.gz --hic-file hic.paired.fq.gz genome.fa
```

Run `python3 scripts/run_pipeline.py -h` for the full option list.

### Pipeline stages

```mermaid
%%{init: {'flowchart': {'rankSpacing': 25, 'nodeSpacing': 15}}}%%
flowchart LR
    G[genome.fa] --> S1[seq_index]
    S1 --> S2[self_align]
    S2 --> S6[hapscotch]
    S1 -. Hi-C only .-> S3[hic_align]
    S3 -.-> S4[hic_convert]
    S4 -.-> S5[contig_ec]
    S5 -.-> S6
    S6 -. Hi-C only .-> S7[yahs_scaffold]
    S6 --> S8[collect_results]
    S7 -.-> S8
    S8 --> R[4.results]
```

| # | Stage | Tool | Produces |
|---|---|---|---|
| 1 | `seq_index` | `seqtools idx` | Genome sequence index (`.idx`) |
| 2 | `self_align` | `selfaln.py` | Self-alignment PAF |
| 3 | `hic_align` | `hicaln.py` | Per-`--hic-file` alignment |
| 4 | `hic_convert` | `hictools convert` | Merged Hi-C BIN |
| 5 | `contig_ec` | `hapcure` | Contig error-correction AGP |
| 6 | `hapscotch` | `hapscotch` | Haplotype-phased scaffolds |
| 7 | `yahs_scaffold` | `seqtools` / `yahs` / `hictools` / `hicmap.py` | Hi-C-rescaffolded haplotypes |
| 8 | `collect_results` | – | Final files under `4.results/` |

### Final outputs

Below is a summary table for the final outputs that can be found in folder `OUTDIR/4.results/`.

* The file `haps.ctg-ec.agp` is only available when running contig error correct.

* Files described in rows 7-10 are only available when running the YaHS module.

| # | File | Description |
|---|---|---|
|1| `haps.ctg-ec.agp` | Contig error correction result | 
|2| `haps.cnt.txt` | Estimated haplotype number/ploidy |
|3| `haps.grp.agp` | Scaffold groups constructed from synteny |
|4| `haps.grp.txt` | Per-contig haplotype assignment |
|5| `haps-[1-9][0-9]*.ctg.agp` | Contig assignment — each haplotype |
|6| `haps-[1-9][0-9]*.ctg.fa.gz` | Contig assignment — each haplotype |
|7| `haps-all.scf.agp` | Scaffolds combining Hi-C — all haplotypes |
|8| `haps-[1-9][0-9]*.scf.agp` | Scaffolds combining Hi-C — each haplotype |
|9| `haps-[1-9][0-9]*.scf.fa.gz` | Scaffolds combining Hi-C — each haplotype |
|10| `haps-all.scf.hic.png` / `.pdf` | Hi-C contact map for the scaffolded assembly |

These files can be combined with `seqtools seq` to generate various FASTA outputs:

**Scaffolds from all haplotypes** (concatenation of `haps.[1-9][0-9]*.fa`):

```bash
# -a for AGP input, -z for BGZIP output
seqtools seq -a -z genome.fa haps-all.scf.agp >haps-all.scf.fa.gz

# equivalent: a '.gz' output suffix also triggers -z
seqtools seq -a -o haps-all.scf.fa.gz genome.fa haps-all.scf.agp
```

**Contig sequences assigned to a given haplotype:**

```bash
seqtools seq genome.fa <(awk '$5==1' haps.grp.txt) >haps-1.ctg.fa
seqtools seq genome.fa <(awk '$5==3' haps.grp.txt) >haps-3.ctg.fa
```

## Running step by step

Each stage can also be run directly. The options shown below are the minimum needed — run any tool with `-h`/`--help` for the full list.

### 1. Build a sequence index

```bash
seqtools idx -o genome.fa.idx genome.fa
```

Builds a genome index file used by several other tools. The genome file should contain contig sequences from all haplotypes. For example, if you use [hifiasm](https://github.com/chhylp123/hifiasm), you could do something like:

```bash
cat hifiasm-hic.hap*.p_ctg.gfa | awk '/^S/{print ">"$2"\n"$3}' >genome.fa
```

If you plan to use Hi-C phasing, it is highly recommended to also run hifiasm with Hi-C data.

### 2. Self-alignment

```bash
FastGA -v -T8 -paf genome.fa >genome.paf
```

Generates an all-vs-all sequence alignment with FastGA, used as input to `hapscotch`. You can gzip the output PAF to save space.

The `scripts/selfaln.py` script wraps this step for the end-to-end pipeline. Since it can be time-consuming, it is recommended to run it yourself and pass the PAF[.gz] file to the pipeline with `--seq-aln` for better control.

### 3. (optional) Align Hi-C reads

```bash
# index the genome
minibwa index -t4 genome.fa genome

# align hic reads in FASTQ format
minibwa map --hic -t16 genome R1.fq.gz R2.fq.gz | \
    samtools fixmate -mpu - - | \
    samtools sort --write-index -l1 -o hic.srt.bam

# mark duplicates
samtools markdup --write-index -c -@16 hic.srt.bam hic.mkdup.bam

# sort by read name
samtools sort -N -@16 -o hic.bam hic.mkdup.bam
```

This is a minimal pipeline for aligning Hi-C reads with `minibwa`, implemented in `scripts/hicaln.py` (used by `scripts/run_pipeline.py`). There are many other Hi-C mapping pipelines you could use: [Arima Genomics](https://github.com/ArimaGenomics/mapping_pipeline), [Omni-C](https://omni-c.readthedocs.io/en/latest/), [HiC-Pro](https://github.com/nservant/HiC-Pro).

If you have multiple Hi-C libraries, repeat this — except for the genome index step — for each library, giving one BAM file per library.

You can also give `minibwa` a single FASTQ file containing interleaved paired-end Hi-C reads or a `stdin` stream input. For example, if you have a name-sorted BAM/CRAM:

```bash
samtools fasta -F0xB00 -n hic-in.bam | \
    minibwa map --hic -t16 genome - | \
    ...
```

The script `scripts/hicaln.py` accepts all of these input types, and a mixture of them.

This step can also be computationally intensive, so it is highly recommended to run it yourself and pass the resulting alignment file(s) to the pipeline with the repeatable `--hic-aln` option.


### 4. (optional) Convert Hi-C alignment(s) to a binary index

```bash
hictools convert -o hicaln-ctg genome.fa.idx hic.1.bam hic.2.bam
```

Merges/converts Hi-C alignments from multiple libraries into a single binary file, used as input by several other tools such as `hapcure` and `hapscotch`. Other input types are also accepted, including BED, PA5, or additional BIN files.

### 5. (optional) Contig error-correction from Hi-C signal

```bash
hapcure -o ctg genome.fa.idx hicaln-ctg.bin
```

Produces `ctg.ec.agp`, which can be passed to `hapscotch` via `-a` option. By default, `hapscotch` runs error-correction automatically when `-a` is not given and `-c` (Hi-C data) is, unless you skip it explicitly with `-E`.

### 6. Haplotype binning/phasing

```bash
hapscotch -o haps -a ctg.ec.agp -c hicaln-ctg.bin -Y genome.fa.idx genome.paf
```

Bins sequences into haplotype groups and builds longer scaffolds using sequence overlaps. With `-c` (Hi-C reads), it also phases haplotypes to minimise inter-haplotype Hi-C signal.

The `-Y` option additionally writes `haps.grp-hic.bin`, used as input for YaHS scaffolding in the next step.

### 7. (optional) Hi-C rescaffolding with YaHS

```bash
seqtools seq -o haps.bbseq.fa -a genome.fa haps.bbseq.agp
## YaHS needs .fai index
seqtools idx -o haps.bbseq.fa.fai haps.bbseq.fa
yahs -a haps.bbscf.agp --no-contig-ec --no-scaffold-ec -o haps haps.bbseq.fa haps.bbseq-hic.bin
```

This builds a backbone sequence representing each scaffold group (which contains sequences from all haplotypes). Hi-C data mapped within each group is compressed and transformed into backbone coordinates. YaHS is finally used to scaffold the backbones into larger scaffolds.

The YaHS AGP output (`haps_scaffolds_final.agp`) and the HapScotch backbone file (`haps.bbpos.txt`) can be combined with `seqtools` to produce various outputs:

**Scaffolds from all haplotypes:**

```bash
# generate an AGP file
seqtools hap -o haps-all.scf.agp haps_scaffolds_final.agp haps.bbpos.txt
# generate a FASTA file using the AGP file
seqtools seq -a -o haps-all.scf.fa.gz genome.fa haps-all.scf.agp
```

**Scaffolds for a given haplotype:**

```bash
seqtools hap -p 2 -o haps-2.scf.agp haps_scaffolds_final.agp haps.bbpos.txt
seqtools seq -a genome.fa haps-2.scf.agp >haps-2.scf.fa
```

**Diagnostic Hi-C contact plots (optional):**

```bash
hictools prepare -a haps-all.scf.agp -n 2000 -o haps-all.scf.hic.txt hicaln-ctg.bin genome.fa.idx
python3 scripts/hicmap.py --png haps-all.scf.hic.png --pdf haps-all.scf.hic.pdf haps-all.scf.hic.txt
```

## Other tools

**`hapcount`** is a standalone command-line tool for ploidy estimation, and is also used internally by `hapscotch`:

```bash
hapcount -o out genome.fa.idx genome.paf
```
