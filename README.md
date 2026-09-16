# HapScotch - haplotype resolved polyploid genome scaffolding

## Overview
HapScotch phases and scaffolds haplotypes in polyploid/multi-haplotype genome
assemblies from all-vs-all sequence self-alignments, optionally refined with
Hi-C data (contig error-correction and Hi-C-guided rescaffolding via
[YaHS](https://github.com/c-zhou/yahs)). The repo builds five C programs -
`hapscotch` (core phasing/scaffolding), `hapcure` (Hi-C contig error-correction),
`hictools` (Hi-C alignment file processing), `seqtools` (sequence/AGP
utilities) and `hapcount` (standalone ploidy estimator) - plus a
`scripts/run_pipeline.py` driver that runs them end to end.

## Installation

### From source
You need a C compiler, GNU make, cmake and zlib/pthread development files.

    git clone --recurse-submodules https://github.com/c-zhou/hapscotch.git
    cd hapscotch
    make
    make install  ## optionally install all executables in ~/bin/

If you already cloned without `--recurse-submodules`, run
`git submodule update --init --recursive` first - this fetches the bundled
[HiGHS](https://github.com/ERGO-Code/HiGHS) ILP solver used for Hi-C phasing.
This produces `hapscotch`, `hapcure`, `hictools`, `seqtools` and `hapcount`.
This still leaves the external pipeline tools (FastGA/minibwa/samtools/YaHS)
for you to install separately - see the next options if you want everything
in one step.

### Conda

* [recipes/conda/](recipes/conda/) is a `conda-build` recipe that
builds hapscotch itself as a conda package (with FastGA/minibwa/samtools/YaHS
pulled in as run dependencies), so a single `conda install` gives you everything:

    ```
    ## create a conda environment with all tools
    conda build recipes/conda/
    conda create -n hapscotch --use-local hapscotch
    ## activate the conda environment
    conda activate hapscotch
    ## run pipeline 
    run_pipeline.py --help
    ```

* [recipes/environment.yml](recipes/environment.yml) sets up a conda environment
with the build toolchain plus FastGA/minibwa/samtools/YaHS and the Python
packages used by `hicmap.py`, then you still build hapscotch itself with `make`:

    ```
    ## create a conda environment with all dependencies
    conda env create -f recipes/environment.yml
    conda activate hapscotch
    ## build hapscotch
    make
    make install  ## optionally install all executables in ~/bin/
    ## run pipeline
    python3 scripts/run_pipeline.py --help
    ```

### Docker / Singularity / Apptainer

* [recipes/Dockerfile](recipes/Dockerfile) builds an image with the full
pipeline (all five binaries, `scripts/*.py`, and FastGA/minibwa/samtools/YaHS):

    ```
    docker build -t hapscotch -f recipes/Dockerfile .
    docker run --rm -v "$PWD":/data -w /data hapscotch run_pipeline.py --help
    ```

* [recipes/Singularity.def](recipes/Singularity.def) is the equivalent
Apptainer/Singularity definition, built the same way from the same
conda-forge/bioconda base:

    ```
    ## or replace apptainer with singularity
    apptainer build hapscotch.sif recipes/Singularity.def
    apptainer exec hapscotch.sif run_pipeline.py --help
    ```

## Running the pipeline
The command `python3 scripts/run_pipeline.py` runs whole pipeline and is resumable:
* `--resume` to resume from the failed step
* `--force-from STEP` to force a rerun from a given STEP
* `--force` to force a rerun from the scratch

Prerequisites to run it:
* The five binaries above, either on `PATH` or pointed to with
  `--seqtools-bin`/`--hictools-bin`/`--hapscotch-bin`/`--hapcure-bin`
* [YaHS](https://github.com/c-zhou/yahs) (`--yahs-bin`), only if you want the
  Hi-C rescaffolding branch
* [FastGA](https://github.com/thegenemyers/FASTGA) (`--fastga-bin`) if you do 
  not already have self-alignments bewteen sequences
* [minibwa](https://github.com/lh3/minibwa) and [samtools](https://github.com/samtools/samtools)
  if you do not already have Hi-C alignment

Basic usage (haploid/no Hi-C):
    
    python3 scripts/run_pipeline.py -o OUTDIR genome.fa

With Hi-C data (enables contig error-correction and, by default, the YaHS
rescaffolding branch):

    python3 scripts/run_pipeline.py -o OUTDIR --hic-aln hic.bam genome.fa
    python3 scripts/run_pipeline.py -o OUTDIR --hic-file r1.fq.gz,r2.fq.gz --hic-file hic.paired.fq.gz genome.fa

Run `python3 scripts/run_pipeline.py -h` for the full option list.

Pipeline stages:

| # | stage           | tool                                     | produces                       |
|---|-----------------|------------------------------------------|--------------------------------|
| 1 | seq_index       | `seqtools idx`                           | genome sequence index (`.idx`) |
| 2 | self_align      | `selfaln.py`                             | self-alignment PAF             |
| 3 | hic_align       | `hicaln.py`                              | per-`--hic-file` alignment     |
| 4 | hic_convert     | `hictools convert`                       | merged Hi-C BIN                |
| 5 | contig_ec       | `hapcure`                                | contig error-correction AGP    |
| 6 | hapscotch       | `hapscotch`                              | haplotype phased scaffolds     |
| 7 | yahs_scaffold   | `seqtools`/`yahs`/`hictools`/`hicmap.py` | Hi-C-rescaffolded haplotypes   |
| 8 | collect_results | -                                        | final files under `4.results`  |

Final outputs (`OUTDIR/4.results`):

    haps.grp.txt               per-contig haplotype assignment
    haps.grp.agp               scaffold groups constructed from synteny
    haps.cnt.txt               estimated haplotype number/ploidy
                               /** addtional files with YaHS module **/
    haps.all.agp               scaffolds combining Hi-C [all haplotypes]
    haps.[1-9][0-9]*.agp       scaffolds combining Hi-C [each haplotype]
    haps.[1-9][0-9]*.fa.gz     scaffolds combining Hi-C [each haplotype]
    haps.all.hic.png/.pdf      Hi-C contact map for the scaffolded assembly

Those files in the result folder can be used combining with `seqtools seq` to generate various types of FASTA output.

* *Generate a FASTA file with scaffolds from all haplotypes, i.e., concatenation of `haps.[1-9][0-9]*.fa`*

    ```
    seqtools seq -a genome.fa haps.all.agp >haps.all.fa   ##need -a option for AGP input
    ```

* *Generate a FASTA file with all sequences assigned to a given haplotype*

    ```
    seqtools seq genome.fa <(awk '{if($3==1) print $1}') >seqs.h1.fa
    seqtools seq genome.fa <(awk '{if($3==3) print $1}') >seqs.h3.fa
    ```

## Running step by step
Each stage can also be run directly. Options shown are the minimum needed;
run any tool with `-h`/`--help` for the full list.

**1. Build a sequence index**

    seqtools idx -o genome.fa.idx genome.fa

**2. Self-alignment** (produces the PAF used by `hapscotch`)

    FastGA -v -T8 -paf genome.fa > genome.paf

**3. (optional) Align Hi-C reads** (produces name-sorted BAM)

**4. (optional) Convert Hi-C alignment(s) to a binary index** (accepts multiple files, merged into one)

    hictools convert -o hicaln genome.fa.idx hic.1.bam hic.2.bam

**5. (optional) Contig error-correction from Hi-C signal**

    hapcure -o ctg genome.fa.idx hicaln.bin

produces `ctg.ec.agp`, fed back into `hapscotch` with `-a`.

**6. Haplotype phasing**

    hapscotch -o haps -a ctg.ec.agp -c hicaln.bin --file-type BIN -Y genome.fa.idx genome.paf

(`-a`/`-c` are optional; `-Y` additionally writes the inputs needed for step 7)

**7. (optional) Hi-C rescaffolding with YaHS**

    seqtools seq -o haps.bbseq.fa -a genome.fa haps.bbseq.agp
    seqtools idx -o haps.bbseq.fa.idx haps.bbseq.fa
    yahs -a haps.bbscf.agp --no-contig-ec --no-scaffold-ec -o haps.bbseq haps.bbseq.fa haps.bbscf-hic.bin
    seqtools hap -o allhaps.scf.agp haps.bbseq_scaffolds_final.agp haps.bbpos.txt
    seqtools seq -o allhaps.scf.fa -a genome.fa allhaps.scf.agp
    ## optionally for diagnostic hic plots
    hictools prepare -a allhaps.scf.agp -n 2000 -o allhaps.scf.hic.txt haps.bbscf-hic.bin genome.fa.idx
    python3 scripts/hicmap.py --png allhaps.scf.hic.png --pdf allhaps.scf.hic.pdf allhaps.scf.hic.txt

## Other tools

***hapcount*** is a standalone command line tool for ploidy-estimation, which is also part of ***hapscotch***

    hapcount -o out genome.fa.idx genome.paf

