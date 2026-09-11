# HapScotch - full pipeline image (hapscotch/hapcure/hictools/seqtools/hapcount
# + FastGA/minibwa/samtools/YaHS + scripts/run_pipeline.py and its helpers).
#
# Build:  docker build -t hapscotch -f recipes/Dockerfile .
# Run:    docker run --rm -v "$PWD":/data -w /data hapscotch run_pipeline.py -h
#         docker run --rm -it -v "$PWD":/data -w /data hapscotch   # interactive shell

FROM condaforge/miniforge3:latest

RUN mamba create -y -n hapscotch -c conda-forge -c bioconda \
        python numpy matplotlib-base make cmake zlib \
        c-compiler cxx-compiler git samtools fastga minibwa yahs \
    && mamba clean -afy
ENV PATH=/opt/conda/envs/hapscotch/bin:$PATH

RUN git clone --recurse-submodules https://github.com/c-zhou/hapscotch.git /opt/hapscotch \
    && cd /opt/hapscotch \
    && make -j"$(nproc)" \
    && mkdir -p /opt/conda/envs/hapscotch/bin/lib \
    && cp hapscotch hapcure hictools seqtools hapcount /opt/conda/envs/hapscotch/bin/ \
    && cp -a lib/. /opt/conda/envs/hapscotch/bin/lib/ \
    && cp scripts/run_pipeline.py scripts/selfaln.py scripts/hicaln.py scripts/hicmap.py \
        /opt/conda/envs/hapscotch/bin/ \
    && chmod +x /opt/conda/envs/hapscotch/bin/*.py \
    && rm -rf /opt/hapscotch

WORKDIR /data
CMD ["/bin/bash"]
