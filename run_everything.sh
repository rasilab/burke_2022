base_folder=$(pwd)
echo "Running analysis in $base_folder folder"
# For singularity only, run this from command line before running this script
# module load Singualarity
# conda activate snakemake 
singularity pull docker://ghcr.io/rasilab/burke_2022:latest

# download SRA annotations and FASTQ files from SRA
cd $base_folder/analysis/barcode_seq/
sh submit_cluster.sh --snakefile download_sra_annotations.smk $@
sh submit_cluster.sh --snakefile download_fastq.smk $@

# run Pool 3 linkage analysis
cd $base_folder/analysis/barcode_seq/pool3_linkage/scripts
sh submit_cluster.sh $@

# run Pool 3 mRNA barcode counting analysis
cd $base_folder/analysis/barcode_seq/pool3_mrna/scripts
sh submit_cluster.sh $@

# run Riboseq analysis 
cd $base_folder/analysis/riboseq/han2020/scripts
sh submit_cluster "--snakefile" make_bowtie2_databases.smk $@
sh submit_cluster "--snakefile" make_motif_lists.smk $@
sh submit_cluster "--snakefile" analyze_riboseq_data.smk $@

