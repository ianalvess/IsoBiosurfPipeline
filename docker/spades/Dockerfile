FROM condaforge/miniforge3:latest

RUN mamba install -y -c bioconda -c conda-forge spades=3.15.5 && \
    mamba clean -afy

WORKDIR /data
ENTRYPOINT ["spades.py"]