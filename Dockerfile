ARG BASEIMAGE=epi2melabs/base-workflow-image:latest
FROM $BASEIMAGE

RUN \
    . $CONDA_DIR/etc/profile.d/mamba.sh \
    && micromamba activate \ 
    && micromamba install -y \
            medaka==1.2.1 \
            minimap2==2.17 \
            racon==1.4.13 \
        -c anaconda -c conda-forge -c bioconda -q -y \
    && fix-permissions $CONDA_DIR \
    && fix-permissions $HOME

USER $WF_UID
WORKDIR $HOME
