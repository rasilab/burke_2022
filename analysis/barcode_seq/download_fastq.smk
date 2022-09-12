import pandas as pd
import itertools as it

sra_annotations = pd.read_table("../../annotations/sra_annotations.tsv")
print(sra_annotations)

container: "../../burke_2022_latest.sif"

rule all:
  input:
    fastq = [f'../../data/fastq/{srr}.fastq' for srr in sra_annotations['srr']]


rule get_fastq:
  """Download fastq from SRA"""
  input: '../../annotations/sra_annotations.tsv'
  output: '../../data/fastq/{srr}.fastq'
  params:
    directory = '../../data/fastq/'
  conda: "sratools"
  threads : 36
  shell:
    """
    set +e # continue if there is an error code
    fasterq-dump --concatenate-reads --outdir {params.directory} --threads {threads} {wildcards.srr}
    exitcode=$?
    if [ $exitcode -eq 1 ]
    then
      exit 1
    else
      exit 0
    fi
    """
