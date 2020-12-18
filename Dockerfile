FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive
LABEL maintainer="Oxford Nanopore Technologies"

ARG WF_USER="epi2melabs"
ARG WF_UID="1000"
ARG WF_GID="100"

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN \
    apt-get update \
    && apt-get install -yq --no-install-recommends \
        wget ca-certificates sudo \ 
    && apt-get clean && rm -rf /var/lib/apt/lists/*

ENV CONDA_DIR=/home/$WF_USER/conda \
    SHELL=/bin/bash \
    WF_USER=$WF_USER \
    WF_UID=$WF_UID \
    WF_GID=$WF_GID
ENV HOME=/home/$WF_USER \
    PATH=$CONDA_DIR/bin:$PATH

COPY fix-permissions /usr/local/bin/fix-permissions
RUN chmod a+rx /usr/local/bin/fix-permissions

RUN \
    echo "auth requisite pam_deny.so" >> /etc/pam.d/su \
    && sed -i.bak -e 's/^%admin/#%admin/' /etc/sudoers \
    && sed -i.bak -e 's/^%sudo/#%sudo/' /etc/sudoers \
    && useradd -m -s /bin/bash -N -u $WF_UID $WF_USER \
    && chmod g+w /etc/passwd \
    && fix-permissions $HOME

USER $WF_UID
WORKDIR $HOME
RUN \
    wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba \
    && ./bin/micromamba shell --help \
    && ./bin/micromamba shell init -s bash -p $CONDA_DIR \
    && source ~/.bashrc \
    && chown $WF_USER:$WF_GID $CONDA_DIR \
    && fix-permissions $CONDA_DIR \
    && fix-permissions $HOME

ENV MAMBA_EXE="$HOME/bin/micromamba"
ENV MAMBA_ROOT_PREFIX="$HOME/conda"

RUN \
    . $CONDA_DIR/etc/profile.d/mamba.sh \
    && micromamba activate \ 
    && micromamba install -y \
            bedtools \
            csvtk \
            curl \
            gffcompare \
            gffread \
            git \
            matplotlib-base \
            minimap2 \
            pandas \
            pychopper \
            requests \
            samtools \
            seaborn-base \
            seqkit \
            snakemake \
            stringtie \
        -c conda-forge -c bioconda -q -y \
    && fix-permissions $CONDA_DIR \
    && fix-permissions $HOME

RUN \
    git clone --depth 1 --single-branch --branch v1.0.0 \
        https://github.com/nanoporetech/pipeline-nanopore-ref-isoforms.git \
    && rm -rf pipeline-nanopore-ref-isoforms/evaluation

COPY run_workflow $CONDA_DIR/bin/
ENTRYPOINT ["/home/epi2melabs/conda/bin/run_workflow"]
CMD ["--help"]

USER $WF_UID
WORKDIR $HOME
