base_folder=$(pwd)
echo "Running analysis in $base_folder folder"
# For singularity only, run this from command line before running this script
# module load Singualarity
# conda activate snakemake 
echo "Pulling singularity container from Github"
singularity pull docker://ghcr.io/rasilab/burke_2022:latest

echo "Downloading SRA annotations from GEO"
cd $base_folder/analysis/barcode_seq/
sh submit_cluster.sh "--snakefile" download_sra_annotations.smk "--forceall" $@
echo "Downloading FASTQ annotations from SRA"
sh submit_cluster.sh "--snakefile" download_fastq.smk $@

echo "Running Pool 3 linkage analysis"
cd $base_folder/analysis/barcode_seq/pool3_linkage/scripts
sh submit_cluster.sh "--snakefile" get_pool3_linkage.smk $@

echo "Running Pool 3 mRNA analysis"
cd $base_folder/analysis/barcode_seq/pool3_mrna/scripts
sh submit_cluster.sh "--snakefile" run_pool3_mrna_barcode_count.smk $@

echo "Running Riboseq analysis"
cd $base_folder/analysis/riboseq/han2020/scripts
echo "	Making Bowtie2 databases"
sh submit_cluster.sh "--snakefile" make_bowtie2_databases.smk $@
echo "	Making motif lists"
sh submit_cluster.sh "--snakefile" make_motif_lists.smk $@
echo "	Processing Riboseq reads"
sh submit_cluster.sh "--snakefile" download_sra_annotations.smk "--forceall" # $@
# sh submit_cluster.sh "--snakefile" analyze_riboseq_data.smk $@