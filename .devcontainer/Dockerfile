# base conda container
FROM continuumio/miniconda3:4.12.0

# Install mamba installer for quick conda installations
RUN conda install -y mamba -c conda-forge

# Install standard python packages in base conda environment
RUN mamba install -y -c conda-forge jupyter ipykernel

# install python packages
RUN mamba install -y -c conda-forge pandas biopython plotnine

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
  bioconductor-rtracklayer

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

# sra tools
RUN mamba create -n sratools
RUN mamba install -y -n sratools -c bioconda -c conda-forge sra-tools=2.11.0

# bowtie2 and samtools
RUN mamba create -n bowtie2_samtools
RUN mamba install -y -n bowtie2_samtools -c bioconda -c conda-forge bowtie2=2.4.5
RUN mamba install -y -n bowtie2_samtools -c bioconda -c conda-forge samtools=1.15.1

# for reading fasta sequence in python
RUN mamba install -y -c bioconda pyfaidx

# genome file
RUN mamba install -n R -y -c bioconda bioconductor-bsgenome.hsapiens.ucsc.hg38

# for reading GFF files in python
RUN conda run pip install gffpandas

# for random resampling of data
RUN mamba install -n R -y -c conda-forge r-rsample

# lab R templates
RUN mamba install -n R -y -c conda-forge r-devtools
RUN conda run -n R Rscript -e "devtools::install_github('rasilab/rasilabRtemplates', ref='7d428f9')"

# for bootstrapping R functions
Run mamba install -n R -y -c conda-forge r-boot
