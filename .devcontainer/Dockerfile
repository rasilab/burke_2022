FROM qmcgaw/latexdevcontainer:latest

USER root
WORKDIR /tmp
# Install pandoc
RUN wget https://github.com/jgm/pandoc/releases/download/2.17/pandoc-2.17-1-amd64.deb
RUN dpkg -i pandoc-2.17-1-amd64.deb

# Install conda and python
# Install conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh -O miniconda.sh
RUN bash miniconda.sh -b -p /opt/conda
# Put conda in PATH 
ENV CONDA_DIR /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
# Set up shell for conda activation
RUN eval "$(conda shell.bash hook)"

# Install pandoc filter extension
RUN conda run -n base pip install pandocfilters

# Update Latex installer
RUN tlmgr repository add ftp://ftp.math.utah.edu/pub/tex/historic/systems/texlive/2021/tlnet-final/
RUN tlmgr option repository ftp://ftp.math.utah.edu/pub/tex/historic/systems/texlive/2021/tlnet-final/
RUN tlmgr update --self

# Install Xelatex
RUN apt install -y libfontconfig1 
# to prevent error messages during xetex install
RUN tlmgr update texlive-scripts
RUN tlmgr install xetex

# Latex packages
RUN tlmgr install setspace \
                  unicode-math \
                  xcolor \
                  booktabs \
                  etoolbox \
                  lineno \
                  caption \
                  float \
                  enumitem

# copy Helvetica font
COPY .install/fonts/ /usr/share/fonts/truetype/helvetica/

# Install mamba installer for quick conda installations
RUN conda install mamba -c conda-forge

# Install Python conda environment
RUN mamba create -y -n py
RUN mamba install -y -n py -c conda-forge -c bioconda \
    jupyter \
    pandas \
    matplotlib \
    ipykernel \
    biopython \
    pysam \
    plotnine

RUN mamba create -y -n snakemake
RUN mamba install -y -n snakemake -c conda-forge -c bioconda snakemake-minimal

# Make snakemake visible to python environment
ENV PATH="$PATH:/opt/conda/envs/snakemake/bin"

# Install R conda environment
RUN mamba create -y -n R

RUN mamba install -y -n R -c conda-forge \
    r-tidyverse \
    r-janitor \
    r-irkernel \
    r-glue \
    r-devtools \
    r-plotrix

RUN mamba install -y -n R -c bioconda -c conda-forge \
  bioconductor-plyranges \
  bioconductor-genomicfeatures \
  bioconductor-rtracklayer \
  bioconductor-flowcore

# Set up R jupyter kernel and make it visible to python
ENV PATH="/opt/conda/envs/py/bin:$PATH"
RUN /opt/conda/envs/R/bin/R -s -e "IRkernel::installspec(sys_prefix = T)"
# Set jupyter data dir for discovering kernels
ENV JUPYTER_DATA_DIR="/opt/conda/envs/py/share/jupyter"
ENV JUPYTER_CONFIG_DIR="/opt/conda/envs/py/share/jupyter"

# install rasilab R templates
RUN /opt/conda/envs/R/bin/R -s -e "devtools::install_github('rasilab/rasilabRtemplates')"

# Make R visible to python environment
ENV PATH="$PATH:/opt/conda/envs/R/bin"

# Make 'py' as default conda environment
RUN sed -i 's/conda activate base/conda activate py/' /root/.bashrc

# List available Jupyter kernels and their locations for use in singularity
RUN jupyter kernelspec list > /jupyter.kernels

# Install inkscape
RUN apt update --allow-releaseinfo-change -y && apt install -y inkscape

RUN mamba install -n py -c conda-forge regex

RUN mamba install -n R -c conda-forge r-ggrepel
RUN mamba install -n R -c conda-forge r-ggthemes
