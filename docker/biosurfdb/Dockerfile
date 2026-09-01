FROM condaforge/miniforge3:latest

RUN conda install -y -n base -c conda-forge mamba && \
    conda clean -afy

RUN mamba create -y -n biosurfdb -c bioconda -c conda-forge prodigal=2.6.3 diamond=2.1.9 && \
    mamba clean -afy

ENV PATH=/opt/conda/envs/biosurfdb/bin:$PATH

WORKDIR /data

ENTRYPOINT ["/bin/bash"]