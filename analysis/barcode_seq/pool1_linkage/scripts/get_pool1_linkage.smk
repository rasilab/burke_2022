"""Workflow for identifying barcodes from linkage sequencing of inserts and barcodes

  :Author: Arvind Rasi Subramaniam
  :Date: 26 Mar 2021
"""

# useful libraries
import os
import pandas as pd
import re


# configuration specific to this analysis
sra_annotations = pd.read_table("../../../../annotations/sra_annotations.tsv")
sample_annotations = pd.read_table("../annotations/sample_annotations.csv", 
                                   sep=",", comment = "#", dtype=object)
fastq_dir = "../../../../data/fastq"                                   
fastq_files = os.listdir(fastq_dir)

container: "../../../../burke_2022_latest.sif"
rule all:
  """List of all files we want at the end
  """
  input:
    [f'../data/bowtie2_reference/{sample_name}.1.bt2' for sample_name in sample_annotations.loc[:, "sample_name"]],
    [f'../data/fastq/{sample_name}.fastq' for sample_name in sample_annotations.loc[:, "sample_name"]], 
    [f'../data/alignments/{sample_name}.bam' for sample_name in sample_annotations.loc[:, "sample_name"]],
    [f'../data/filtered_alignments/{sample_name}.tsv.gz' for sample_name in sample_annotations.loc[:, "sample_name"]],
    [f'../data/insert_barcode_counts/{sample_name}.tsv.gz' for sample_name in sample_annotations.loc[:, "sample_name"]],
    [f'../data/ref_vs_ref_alignments/{sample_name}/alignment_barcode1.bam' for sample_name in sample_annotations.loc[:, "sample_name"]],
    [f'../data/filtered_barcodes/{sample_name}.tsv.gz' for sample_name in sample_annotations.loc[:, "sample_name"]],


def get_fastq_file_for_sample_name(wildcards):
  """This function gets the R1 or R2 file depending on the `insert_read` column of `sample_annotations`
  """
  sample_id = sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'sample_id'].item()
  srr = sra_annotations.loc[sra_annotations['sample_id'] == sample_id, 'srr'].item()
  filename = f'{fastq_dir}/{srr}.fastq'
  return filename


rule generate_reference_for_alignment:
  """Generate fasta file from library annotations file for bowtie index reference creation
  """
  input:
    '../annotations/insert_annotations.tsv',
  output:
    '../data/bowtie2_reference/{sample_name}.fasta',
    '../annotations/insert_annotations/{sample_name}.tsv'
  log:
    '../data/bowtie2_reference/{sample_name}_reference_creation.log'
  container: "docker://ghcr.io/rasilab/r:1.0.0"
  shell:
    """
    mkdir -p ../annotation/insert_annotations
    Rscript create_reference_for_alignment.R {input} {output} &> {log}
    """

rule fastq_rename:
  """Rename header of reads in SRA fastq files to consecutive integers for easier comparison
  """
  input:
    fastq = get_fastq_file_for_sample_name
  output: '../data/fastq/{sample_name}.fastq'
  log:
    '../data/fastq/{sample_name}.log'
  container: "docker://ghcr.io/rasilab/python:1.0.0"
  shell:
    """
    sed 's/SRR[[:digit:]]\+.\([[:digit:]]\+\).*/\\1/g' {input.fastq} > {output}
    """


rule generate_bowtie_reference:
  """Generate bowtie2 reference from fasta reference
  """
  input:
    '../data/bowtie2_reference/{sample_name}.fasta'
  output:
    '../data/bowtie2_reference/{sample_name}.1.bt2'
  log:
    '../data/bowtie2_reference/{sample_name}.bowtie2build.log'
  container: "docker://ghcr.io/rasilab/bowtie2:2.4.5"
  shell:
    'bowtie2-build {input} ../data/bowtie2_reference/{wildcards.sample_name} &> {log}'


rule align:
  """Align R1 or R2 reads against reference of library sequences
  """
  input: '../data/fastq/{sample_name}.fastq'
  output: '../data/alignments/{sample_name}.sam'
  log:
    '../data/alignments/{sample_name}.bowtie2align.log'
  threads: 36
  params:
    trim5 = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'trim5'].item(),
    trim3 = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'trim3'].item(),
    reference = '../data/bowtie2_reference/{sample_name}'
  container: "docker://ghcr.io/rasilab/bowtie2:2.4.5"
  shell:
    """
    bowtie2 \
    -x {params.reference} \
    -U {input} \
    -N 1 \
    -L 22 \
    --end-to-end \
    --trim5 {params.trim5} \
    --trim3 {params.trim3} \
    --threads {threads} \
    --nofw \
    --no-unal \
    1> {output} \
    2> {log}
    """


rule sort_and_index_bam:
  """Convert SAM alignments to sorted and indexed BAM alignments
  """
  input:
    '../data/alignments/{sample_name}.sam'
  output:
    bam = '../data/alignments/{sample_name}.bam',
    bam_index = '../data/alignments/{sample_name}.bam.bai'
  params:
    unsorted_bam = '../data/alignments/{sample_name}.unsorted.bam'
  threads: 36
  log:
    '../data/alignments/{sample_name}.samtools.log'
  container: "docker://ghcr.io/rasilab/bowtie2:2.4.5"
  shell:
    """
    # convert to BAM
    samtools view -@ {threads} -b {input} > {params.unsorted_bam} 2> {log}
    # sort
    samtools sort -@ {threads} {params.unsorted_bam} > {output.bam} 2>> {log}
    # index
    samtools index -@ {threads} {output.bam} 2>> {log}
    """


rule filter_alignments:
  """Filter alignments to remove duplicates and set maximum mismatches
  """
  input:
    bam = '../data/alignments/{sample_name}.bam',
    Rscript = 'filter_alignments.R'
  output:
    tsv = '../data/filtered_alignments/{sample_name}.tsv.gz',
  log:
    '../data/filtered_alignments/{sample_name}.log',
  container: "docker://ghcr.io/rasilab/r:1.0.0"
  shell:
    'Rscript {input.Rscript} {input.bam} {output.tsv}'


rule count_insert_barcode_pairs:
  """Parse barcodes for each filtered alignment from R1 fastq file and tabulate count for each insert-barcode
  """
  input:
    filtered_alignments = '../data/filtered_alignments/{sample_name}.tsv.gz',
    fastq = '../data/fastq/{sample_name}.fastq'
  output:
    tsv = '../data/insert_barcode_counts/{sample_name}.tsv.gz',
  params:
    barcode1_read = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'barcode1_read'].item(),
    barcode1_start = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'barcode1_start'].item(),
    barcode1_length = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'barcode1_length'].item(),
  log:
    '../data/insert_barcode_counts/{sample_name}.log',
  container: "docker://ghcr.io/rasilab/python:1.0.0"
  shell:
    """
    python count_barcode_insert_pairs.py {input.fastq} {input.filtered_alignments} {output.tsv} \
        {params.barcode1_read} {params.barcode1_start} {params.barcode1_length} \
        &> {log}
    """


rule align_barcodes_1_against_themselves:
  """Align barcodes 1 against themselves to find multialigners
  """
  input:
    '../data/insert_barcode_counts/{sample_name}.tsv.gz'
  output:
    sam = temp('../data/ref_vs_ref_alignments/{sample_name}/alignment_barcode1.sam'),
    bam = '../data/ref_vs_ref_alignments/{sample_name}/alignment_barcode1.bam',
    fasta = '../data/ref_vs_ref_alignments/{sample_name}/reference_barcode1.fasta',
  log:
    align = '../data/ref_vs_ref_alignments/{sample_name}/align_barcode1.log',
    build = '../data/ref_vs_ref_alignments/{sample_name}/build_barcode1.log',
  params:
    bowtie_index = '../data/ref_vs_ref_alignments/{sample_name}/reference_barcode1'
  threads: 36
  container: "docker://ghcr.io/rasilab/bowtie2:2.4.5"
  shell:
    """
    # write the input file to a fasta file of barcodes with name as barcode_num col from input_file
    zcat {input} | awk 'NR > 1 {{print ">" $1 "\\n" $2}}' > {output.fasta}
    # create a bowtie reference of the barcodes
    bowtie2-build {output.fasta} {params.bowtie_index} 2> {log.build}
    # align against itself
    bowtie2 --threads {threads} -L 19 -N 1 --all --norc --no-unal -f -x {params.bowtie_index} -U {output.fasta} > {output.sam}  2> {log.align}
    # convert to BAM
    samtools view -@ {threads} -b {output.sam} > {output.bam}.tmp
    # sort
    samtools sort -@ {threads} {output.bam}.tmp > {output.bam}
    sleep 10
    # index
    samtools index -@ {threads} {output.bam}
    # remove unsorted bam
    rm {output.bam}.tmp
    """


rule filter_barcodes:
  """Filter barcodes to remove clashes and sequencing errors and produce a final list
  """
  input:
    bam1 = '../data/ref_vs_ref_alignments/{sample_name}/alignment_barcode1.bam',
    counts = '../data/insert_barcode_counts/{sample_name}.tsv.gz',
    Rscript = 'filter_barcodes.R',
  output:
    '../data/filtered_barcodes/{sample_name}.tsv.gz',
  container: "docker://ghcr.io/rasilab/r:1.0.0"
  shell:
    'Rscript {input.Rscript} {input.bam1} {input.counts} {output}'
