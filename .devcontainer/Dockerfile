# base conda container
FROM continuumio/miniconda3:4.12.0

# Install mamba installer for quick conda installations
RUN conda install -y mamba -c conda-forge

# Install standard python packages in base conda environment
RUN mamba install -y -c conda-forge jupyter ipykernel

# Create a separate R conda environment 
RUN mamba create -y -n R
# tidyverse and R jupyter kernel
RUN mamba install -y -n R -c conda-forge \
    r-tidyverse \
    r-janitor \
    r-irkernel \
    r-plotrix

# bioconductor packages
RUN mamba install -y -n R -c bioconda -c conda-forge \
  bioconductor-plyranges \
  bioconductor-genomicfeatures \
  bioconductor-rtracklayer \
  bioconductor-bsgenome.hsapiens.ucsc.hg38    

# Make R kernel visible to Jupyter
RUN /opt/conda/envs/R/bin/R -s -e "IRkernel::installspec(sys_prefix = T)"

# packages for downloading SRA annotations
RUN mamba install -y -n R -c bioconda -c conda-forge \
  bioconductor-geoquery \
  r-reutils \
  r-xml2 \
  r-r.utils \
  r-glue

# workflow manager and dataframe package
RUN mamba install -y -n base -c bioconda -c conda-forge snakemake-minimal pandas

RUN mamba install -y -c conda-forge singularity
