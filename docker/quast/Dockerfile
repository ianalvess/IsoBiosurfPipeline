FROM condaforge/miniforge3:latest

RUN conda install -y -n base -c conda-forge mamba && \
    conda clean -afy

RUN mamba create -y -n quast -c bioconda -c conda-forge quast=5.2.0 && \
    mamba clean -afy

ENV PATH=/opt/conda/envs/quast/bin:$PATH

WORKDIR /output

ENTRYPOINT ["quast.py"]