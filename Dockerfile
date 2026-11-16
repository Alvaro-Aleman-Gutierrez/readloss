FROM mambaorg/micromamba:1.5.8
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/readloss
WORKDIR /opt/readloss
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 scipy pytest diamond=2.1.10 "sra-tools>=3.0" "blast>=2.12" matplotlib-base "snakemake-minimal>=8" pip && \
    micromamba run -n base pip install --no-deps . && \
    micromamba clean -a -y
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "readloss"]
