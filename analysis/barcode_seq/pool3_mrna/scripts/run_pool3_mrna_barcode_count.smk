"""Workflow for identifying counting barcodes in dual barcode mRNA sequencing

  :Author: Arvind Rasi Subramaniam
  :Date: 27 Mar 2021
"""

# useful libraries
import os
import pandas as pd
import re
import itertools as it


# configuration specific to this analysis
insert_barcode_linkage_folder = '../../pool3_linkage/data/filtered_barcodes'
sra_annotations = pd.read_table("../../../../annotations/sra_annotations.tsv")
sample_annotations = pd.read_table("../annotations/sample_annotations.csv", 
                                   sep=",", comment = "#", dtype=object)
fastq_dir = "../../../../data/fastq"                                   
fastq_files = os.listdir(fastq_dir)


def get_fastq_file_for_sample_name(wildcards):
  """This function gets the R1 or R2 file depending on the `insert_read` column of `sample_annotations`
  """
  sample_id = sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'sample_id'].item()
  srr = sra_annotations.loc[sra_annotations['sample_id'] == sample_id, 'srr'].item()
  filename = f'{fastq_dir}/{srr}.fastq'
  return filename


def get_bowtie_reference_for_sample_name(wildcards, **kwargs):
  print(wildcards)
  linkage_ref = sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'linkage_ref'].item()
  try:
      barcode_direction = kwargs['barcode_direction']
  except KeyError:
      print("Need a barcode direction (1 or 2) to get the correct bowtie reference file.")

  index = f'../data/bowtie_references/{linkage_ref}/reference_barcode_{barcode_direction}.1.bt2'
  return index

container: "../../../../burke_2022_latest.sif"

# these rules are run locally
localrules: all

# Rules ----------------------------------------------------------------------

rule all:
  """List of all files we want at the end
  """
  input:
    # summarized_counts = ['../tables/sample_insert_barcode_counts.tsv.gz'],
    barcode_counts = [f'../data/barcode_counts/{sample_name}_barcode_{barcode_dir}.tsv.gz' 
                      for (sample_name, barcode_dir) in it.product(sample_annotations['sample_name'].tolist(), [1])],
    barcode_alignments = [f'../data/alignments/{sample_name}_barcode_{barcode_dir}.bam'
                      for (sample_name, barcode_dir) in it.product(sample_annotations['sample_name'].tolist(), [1])],
    bowtie_references = ['../data/bowtie_references/' + file_name.split(".")[0] + '/reference_barcode_1.fasta' 
        for file_name in os.listdir(insert_barcode_linkage_folder) if file_name.endswith(".tsv.gz")],


rule create_bowtie_reference:
  """Create bowtie2 reference for aligning dual barcodes from linkage seq
  """
  input:
    insert_barcode_linkage_folder + '/{sample_name}.tsv.gz'
  output:
    fasta1 = '../data/bowtie_references/{sample_name}/reference_barcode_1.fasta',
    bowtie_ref1 = ([('../data/bowtie_references/{sample_name}/reference_barcode_1.' + str(n) + '.bt2') for n in range(1,5)] + 
        [('../data/bowtie_references/{sample_name}/reference_barcode_1.rev.' + str(n) + '.bt2') for n in range(1,3)]),
  params:
    reference1 = '../data/bowtie_references/{sample_name}/reference_barcode_1',
  log:
    reference1 = '../data/bowtie_references/{sample_name}/reference_barcode_1.create_reference.log',
  shell:
    """
    # write the input file to a fasta file of barcodes that are serially numbered
    conda run -n R Rscript create_bowtie_reference.R {input} {output.fasta1}
    # create a bowtie reference of the barcodes
    conda run -n bowtie2_samtools bowtie2-build {output.fasta1} {params.reference1} &> {log.reference1}
    """


def get_reference_name_from_bt2(wildcards, input):
  return input.bowtie_reference.replace('.1.bt2', '')


rule align:
  """Align reads to reference barcode_1
  """
  input:
    fastq = get_fastq_file_for_sample_name,
    bowtie_reference = lambda wildcards: get_bowtie_reference_for_sample_name(wildcards, barcode_direction=1)
  output:
    sam = '../data/alignments/{sample_name}_barcode_1.sam'
  log:
    '../data/alignments/{sample_name}_barcode_1.bowtie2.log'
  threads: 36
  params:
    ref = get_reference_name_from_bt2,
    trim5 = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'barcode1_trim5'].item(),
    trim3 = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'barcode1_trim3'].item(),
    seed_length = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'barcode1_length'].item(),
  conda: "bowtie2_samtools"
  shell:
    """
    bowtie2 \
    -L {params.seed_length} \
    --trim5 {params.trim5} \
    --trim3 {params.trim3} \
    --threads {threads} \
    -N 1 --norc --no-unal \
    -x {params.ref} \
    -U {input.fastq} > {output.sam} 2> {log}
    """


rule sort_and_index_bam:
  """Sort and index alignments
  """
  input:
    sam = '../data/alignments/{sample_name}.sam'
  output:
    bam = '../data/alignments/{sample_name}.bam',
    bai = '../data/alignments/{sample_name}.bam.bai',
  threads: 36
  conda: "bowtie2_samtools"
  shell:
    """
    # convert to BAM
    samtools view -@ {threads} -b {input.sam} > {output.bam}.tmp;
    # sort
    samtools sort -@ {threads} {output.bam}.tmp > {output.bam};
    # index
    samtools index -@ {threads} {output.bam}
    # remove unsorted bam
    rm {output.bam}.tmp
    """


rule count:
  """Count the number of reads aligning to each reporter barcode-R1 barcode combination
  """
  
  # The `idxstats` command counts the number of reads per barcode
  # The `awk` command below adds header, removes zero count entries, and sorts by read counts in reverse.
  input:
    '../data/alignments/{sample_name}.bam'
  output:
    '../data/barcode_counts/{sample_name}.tsv.gz'
  log:
    '../data/barcode_counts/{sample_name}.count.log'
  conda: "bowtie2_samtools"
  threads: 8
  shell:
    """
    export TMPDIR=$(pwd) # necessary for sort command below to run inside singularity
    samtools idxstats {input} 1> {wildcards.sample_name}.tmp 2> {log}
    awk 'BEGIN {{OFS="\\t"; print "barcode","count";}} $3 != "0" {{print $1,$3 | "sort -k2nr -k1,1n"}}' {wildcards.sample_name}.tmp | gzip -c > {output} 2> {log}
    rm {wildcards.sample_name}.tmp
    """
