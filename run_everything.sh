base_folder=$(pwd)
echo "Running analysis in $base_folder folder"
cd scripts/
# For singularity only
# module load Singualarity
# conda activate snakemake 
# singularity pull docker://ghcr.io/rasilab/burke_2022:latest

# to run inside singularity container, prepend 
# JUPYTER_DATA_DIR=$PWD
# JUPYTER_CONFIG_DIR=$PWD
# singularity exec --bind /fh burke_2022_latest.sif followed by the commadn below
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir download_sra_annotations.ipynb

# download FASTQ files from SRA
# sh submit_cluster.sh

# run Pool 3 linkage analysis
cd $base_folder
cd scripts/barcode_seq/pool3_linkage/scripts
sh ../../../submit_cluster.sh

# run Pool 3 mRNA barcode counting analysis
cd $base_folder
cd scripts/barcode_seq/pool3_linkage/scripts
# sh ../../../submit_cluster.sh