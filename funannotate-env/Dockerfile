FROM condaforge/miniforge3:26.7.2-0

# procps: Nextflow's trace/report needs `ps` inside the container
RUN apt-get update && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

COPY funannotate.yml /tmp/funannotate.yml
RUN mamba env create -n funannotate -f /tmp/funannotate.yml \
    && mamba clean -afy && rm /tmp/funannotate.yml

ENV CONDA_PREFIX=/opt/conda/envs/funannotate \
    PATH=/opt/conda/envs/funannotate/bin:$PATH \
    AUGUSTUS_CONFIG_PATH=/opt/conda/envs/funannotate/config \
    AUGUSTUS_SCRIPTS_PATH=/opt/conda/envs/funannotate/bin/ \
    AUGUSTUS_BIN_PATH=/opt/conda/envs/funannotate/bin/ \
    QUARRY_PATH=/opt/conda/envs/funannotate/opt/codingquarry-2.0/QuarryFiles \
    EVM_HOME=/opt/conda/envs/funannotate/opt/evidencemodeler-1.1.1 \
    PASAHOME=/opt/conda/envs/funannotate/opt/pasa-2.5.3 \
    TRINITY_HOME=/opt/conda/envs/funannotate/bin \
    JAVA_HOME=/opt/conda/envs/funannotate/lib/jvm \
    JAVA_LD_LIBRARY_PATH=/opt/conda/envs/funannotate/lib/jvm/lib/server \
    HOME=/tmp \
    MPLCONFIGDIR=/tmp
# FUNANNOTATE_DB is deliberately NOT set: the database is mounted at run time.
