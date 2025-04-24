import pandas as pd
import itertools as it

sra_annotations = pd.read_table("../../annotations/sra_annotations.tsv")[-4:]

rule all:
  input:
    fastq = [f'../../data/fastq/{srr}.fastq' for srr in sra_annotations['srr']]


rule get_fastq:
  """Download fastq from SRA"""
  input: '../../annotations/sra_annotations.tsv'
  output: '../../data/fastq/{srr}.fastq'
  params:
    directory = '../../data/fastq/'
  threads : 36
  container: 'docker://ghcr.io/rasilab/sratools:3.0.8'
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
