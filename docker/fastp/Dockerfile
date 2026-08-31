FROM condaforge/miniforge3:latest

RUN mamba install -y -c bioconda -c conda-forge fastp=0.23.4 && \
    mamba clean -afy

WORKDIR /data
ENTRYPOINT ["fastp"]