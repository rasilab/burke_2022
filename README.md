# A Nascent Peptide Code for Translational Control of mRNA Stability in Human Cells <!-- omit in toc -->

**Phil Burke**<sup>1,2</sup>, **Heungwon Park**<sup>1</sup>, **Arvind Rasi Subramaniam**<sup>1,†</sup>

<sup>1</sup> Basic Sciences Division and Computational Biology Section of the
Public Health Sciences Division, Fred Hutchinson Cancer Center, Seattle
<br/>
<sup>2</sup> Microbiology Graduate Program, University of Washington, Seattle <br/>

<sup>†</sup>Corresponding author: A.R.S: <rasi@fredhutch.org>

Nature Communications 13, 6829 (2022); doi: https://doi.org/10.1038/s41467-022-34664-0

**Contents**

- [Abstract](#abstract)

## Abstract

Stability of eukaryotic mRNAs is associated with their codon, amino acid, and GC content.
Yet, coding sequence motifs that predictably alter mRNA stability in human cells remain poorly defined. 
Here, we develop a massively parallel assay to measure mRNA effects of thousands of synthetic and endogenous coding sequence motifs in human cells. 
We identify several families of simple dipeptide repeats whose translation triggers mRNA destabilization.
Rather than individual amino acids, specific combinations of bulky and positively charged amino acids are critical for the destabilizing effects of dipeptide repeats.
Remarkably, dipeptide sequences that form extended β strands *in silico* and *in vitro* slowdown ribosomes and reduce mRNA levels *in vivo*. 
The resulting nascent peptide code underlies the mRNA effects of hundreds of endogenous peptide sequences in the human proteome. 
Our work suggests an intrinsic role for the ribosome as a selectivity filter against the synthesis of bulky and aggregation-prone peptides.


## Software Installation

- Get the docker image from this repo:

```
docker pull ghcr.io/rasilab/burke_2022:latest
```

- Open this repo inside VScode remote development container. It should automatically pick up the above image based on [](./.devcontainer/devcontainer.json).

- Run the analysis workflow from the command line inside the VScode container:

```
cd scripts/
snakemake -np # dry run
snakemake -p --cores=all --use-conda
```

- To run this on a cluster with singularity containers, do:

```
grabnode # for fred hutch only, get max 36 nodes, max 720G memory, for **1** day, with no GPU
module load singularity # for fred hutch cluster
singularity pull docker://ghcr.io/rasilab/burke_2022:latest
conda activate snakemake # this is a minimal conda env that has snakemake-minimal and pandas for invoking snakefile
# for interactive singularity debugging on cluster, run following command
# singularity exec -B /fh:/fh burke_2022_latest.sif /bin/bash
cd $DIRECTORY_OF_INTEREST
sh submit_local.sh # adjust this script to mount the parent directory
# sh submit_cluster.sh # uses SLURM to run the containers in parallel, adjust script to mount parent directory
```
cd analysis
