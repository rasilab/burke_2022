base_folder=$(pwd)
echo "Running analysis in $base_folder folder"
## For singularity only, run this from command line before running this script
# module load Singualarity
# conda activate snakemake 
# echo "Pulling singularity container from Github"
# singularity pull docker://ghcr.io/rasilab/burke_2022:latest

# echo "Downloading SRA annotations from GEO"
# cd $base_folder/analysis/barcode_seq/
# sh submit_cluster.sh "--snakefile" download_sra_annotations.smk "--forceall" $@
# echo "Downloading FASTQ annotations from SRA"
# sh submit_cluster.sh "--snakefile" download_fastq.smk $@

echo "Running Pool 1 linkage analysis"
cd $base_folder/analysis/barcode_seq/pool1_linkage/scripts
sh submit_local.sh "--snakefile" get_pool1_linkage.smk $@

echo "Running Pool 1 mRNA analysis"
cd $base_folder/analysis/barcode_seq/pool1_mrna/scripts
sh submit_local.sh "--snakefile" run_pool1_mrna_barcode_count.smk $@

echo "Running Pool 1 FACS-seq analysis"
cd $base_folder/analysis/barcode_seq/pool1_facs_seq/scripts
sh submit_local.sh "--snakefile" run_pool1_facs_seq_barcode_count.smk $@

echo "Running Pool 2 linkage analysis"
cd $base_folder/analysis/barcode_seq/pool2_linkage/scripts
sh submit_local.sh "--snakefile" get_pool2_linkage.smk $@

echo "Running Pool 2 mRNA analysis"
cd $base_folder/analysis/barcode_seq/pool2_mrna/scripts
sh submit_local.sh "--snakefile" run_pool2_mrna_barcode_count.smk $@

echo "Running Pool 3 linkage analysis"
cd $base_folder/analysis/barcode_seq/pool3_linkage/scripts
sh submit_local.sh "--snakefile" get_pool3_linkage.smk $@

echo "Running Pool 3 mRNA analysis"
cd $base_folder/analysis/barcode_seq/pool3_mrna/scripts
sh submit_local.sh "--snakefile" run_pool3_mrna_barcode_count.smk $@

# echo "Running Pool 4 linkage analysis"
# cd $base_folder/analysis/barcode_seq/pool4_linkage/scripts
# sh submit_local.sh "--snakefile" get_pool4_linkage.smk $@

# echo "Running Pool 4 mRNA analysis"
# cd $base_folder/analysis/barcode_seq/pool4_mrna/scripts
# sh submit_local.sh "--snakefile" run_pool4_mrna_barcode_count.smk $@

# echo "Running Riboseq analysis"
# cd $base_folder/analysis/riboseq/han2020/scripts
# echo "	Making Bowtie2 databases"
# sh submit_cluster.sh "--snakefile" make_bowtie2_databases.smk $@
# echo "	Making motif lists"
# sh submit_cluster.sh "--snakefile" make_motif_lists.smk $@
# echo "	Processing Riboseq reads"
# sh submit_cluster.sh "--snakefile" download_sra_annotations.smk $@
# sh submit_cluster.sh "--snakefile" analyze_riboseq_data.smk $@