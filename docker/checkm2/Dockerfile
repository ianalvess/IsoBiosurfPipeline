FROM condaforge/miniforge3:latest

RUN conda install -y -n base -c conda-forge mamba && \
    conda clean -afy

RUN mamba create -y -n checkm2 -c bioconda -c conda-forge checkm2 && \
    mamba clean -afy

ENV PATH=/opt/conda/envs/checkm2/bin:$PATH

WORKDIR /data

ENTRYPOINT ["checkm2"]