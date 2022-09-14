"""Workflow for ribosome profiling analysis
"""

container:  '../../../../burke_2022_latest.sif'

import pandas as pd
import tabulate

# configuration specific to this analysis
study_annotations = pd.read_table("../annotations/geo_accession_numbers.csv",
                                  sep=",", comment="#")
merged_annotations = pd.read_table("../annotations/sra_annotations.tsv", sep="\t")

# these rules are run locally
localrules: all

# Rules ----------------------------------------------------------------------

rule all:
    """List of all files we want at the end"""
    input:
        fastq = [f'../data/fastq/{srr}.fastq' for srr in merged_annotations.loc[:, "srr"]],
        transcript_alignments = [f'../data/alignments/{srr}.transcripts.bam' for srr in merged_annotations.loc[:, 'srr']],
        transcript_coverage = [f'../data/coverage/{gsm}.transcripts.bedGraph.gz' for gsm in set(merged_annotations.loc[:, "gsm"])],
        transcript_counts = [f'../data/gene_counts/{gsm}.transcripts.tsv.gz' for gsm in set(merged_annotations.loc[:, "gsm"])],
        codon_density = [f'../data/codon_density/{gsm}.tsv.gz' for gsm in set(merged_annotations.loc[:, "gsm"])],
        vk_type_density = [f'../data/vk_type_density/{gsm}.tsv.gz' for gsm in set(merged_annotations.loc[:, "gsm"])],
        vk_type_plots = "plot_vk_type_density_han2020.nbconvert.ipynb",


rule get_fastq:
    """Download fastq from SRA"""
    input:
    output:
        '../data/fastq/{srr}.fastq'
    params:
        directory = '../data/fastq'
    log:
        '../data/fastq/{srr}.log'
    threads:
        36
    conda: "sratools"
    shell:
        """
        set +e # continue if there is an error code
        fasterq-dump --concatenate-reads --outdir {params.directory} --threads {threads} {wildcards.srr} 2> {log}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
          exit 1
        else
          exit 0
        fi
        """


def get_trim_5(wildcards, output):
    trim5 = merged_annotations.loc[merged_annotations['srr']
                                   == wildcards.srr, 'trim5'].item()
    if pd.isnull(trim5):
        return '0'
    trim_condition = merged_annotations.loc[merged_annotations['srr']
                                            == wildcards.srr, 'trim_condition'].item()
    sample_name = merged_annotations.loc[merged_annotations['srr']
                                         == wildcards.srr, 'sample_name'].item()
    if ((not pd.isnull(trim_condition)) and (re.search(trim_condition, sample_name))):
        return f'{int(trim5)}'
    else:
        return '0'


def get_trim_3(wildcards, output):
    trim3 = merged_annotations.loc[merged_annotations['srr']
                                   == wildcards.srr, 'trim3'].item()
    if pd.isnull(trim3):
        return '0'
    trim_condition = merged_annotations.loc[merged_annotations['srr']
                                            == wildcards.srr, 'trim_condition'].item()
    sample_name = merged_annotations.loc[merged_annotations['srr']
                                         == wildcards.srr, 'sample_name'].item()
    if ((not pd.isnull(trim_condition)) and (re.search(trim_condition, sample_name))):
        return f'{int(trim3)}'
    # trim the total RNA sample in Khajuria 2018 by 22 nt to match mean RPF length
    if (('rna' in merged_annotations.loc[merged_annotations['srr'] == wildcards.srr, 'sample_name'].item()) and ('khajuria' in merged_annotations.loc[merged_annotations['srr'] == wildcards.srr, 'study'].item())):
        return '22'
    else:
        return '0'


rule trim_linker:
    """Remove extra sequences before aligning"""
    input:
        '../data/fastq/{srr}.fastq'
    params:
      adapter = lambda wildcards: merged_annotations.loc[merged_annotations['srr'] == wildcards.srr, 'adapter'].item(),
      trim5 = get_trim_5,
      trim3 = get_trim_3
    output:
        temp('../data/trim/{srr}.fastq')
    log:
        '../data/trim/{srr}.log'
    conda: "cutadapt"
    shell:
        """
        cutadapt \
        --adapter={params.adapter} \
        --cut={params.trim5} \
        --cut=-{params.trim3} \
        --minimum-length=22 \
        --match-read-wildcards \
        --output {output} \
        {input} \
        &> {log}
        """


def get_input_for_rrna_subtraction(wildcards):
    trim = merged_annotations.loc[merged_annotations['srr']
                                  == wildcards.srr, 'trim'].item()
    trim_condition = merged_annotations.loc[merged_annotations['srr']
                                            == wildcards.srr, 'trim_condition'].item()
    sample_name = merged_annotations.loc[merged_annotations['srr']
                                         == wildcards.srr, 'sample_name'].item()
    if trim == 'yes':
        if ((pd.isnull(trim_condition)) or (re.search(trim_condition, sample_name))):
            return f'../data/trim/{wildcards.srr}.fastq'
    return f'../data/fastq/{wildcards.srr}.fastq'


rule remove_rrna:
    """Remove contaminant reads aligning to ribosomal rRNA"""
    input:
        get_input_for_rrna_subtraction
    params:
      rrna_index = lambda wildcards: merged_annotations.loc[merged_annotations['srr'] == wildcards.srr, 'rrna_index'].item()
    output:
        temp('../data/norrna/{srr}.fastq')
    log:
        '../data/norrna/{srr}.log'
    conda: "bowtie2_samtools"
    threads: 36
    shell:
        """
        echo "start" &> {log}
        bowtie2 \
        --threads {threads} \
        --un {output} \
        -x {params.rrna_index} \
        -U {input} \
        1> /dev/null \
        2> {log} 
        echo "stop" &> {log}
        """


rule align_transcripts:
    """Align against canonical CCDS transcripts"""
    input:
        '../data/norrna/{srr}.fastq'
    params:
        transcript_index = lambda wildcards: merged_annotations.loc[merged_annotations['srr'] == wildcards.srr, 'transcript_index'].item(),
        rc_or_fw = lambda wildcards: '--nofw' if (('rna' in merged_annotations.loc[merged_annotations['srr'] == wildcards.srr, 'sample_name'].item(
        )) and ('khajuria' in merged_annotations.loc[merged_annotations['srr'] == wildcards.srr, 'study'].item())) else '--norc'
    output:
        '../data/alignments/{srr}.transcripts.sam'
    log:
        '../data/alignments/{srr}.transcripts.alignment.log'
    conda: "bowtie2_samtools"
    threads: 36
    shell:
        """
        bowtie2 \
        {params.rc_or_fw} \
        --threads {threads} \
        -x {params.transcript_index} \
        -U {input} \
        1> {output} \
        2> {log} 
        """


rule sort_and_index_transcript_alignments:
    """Convert SAM alignments to sorted and indexed BAM alignments for transcripts"""
    input:
        '../data/alignments/{srr}.transcripts.sam'
    output:
        unsorted_bam = temp(
            '../data/alignments/{srr}.unsorted_transcripts.bam'),
        bam = '../data/alignments/{srr}.transcripts.bam',
        bam_index = '../data/alignments/{srr}.transcripts.bam.bai'
    log:
        '../data/alignments/{srr}.transcripts.samtools.log'
    conda: "bowtie2_samtools"
    threads: 36
    shell:
        """
        samtools view -@ {threads} -b {input} > {output.unsorted_bam} 2> {log}
        samtools sort -@ {threads} {output.unsorted_bam} > {output.bam} 2>> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}
        """


rule calculate_transcriptomic_coverage:
    """Sum the P-site counts at each transcript location"""
    params:
        study = lambda wildcards: merged_annotations.loc[merged_annotations['gsm'] == wildcards.gsm, 'study'].tolist()[0],
    input:
        trim_distances = '../annotations/p_site_trim_distances.csv',
        alignments = lambda wildcards: [f'../data/alignments/{srr}.transcripts.bam' for srr in merged_annotations.loc[merged_annotations['gsm'] == wildcards.gsm, 'srr'].tolist()]
    output:
        '../data/coverage/{gsm}.transcripts.bedGraph.gz'
    log:
        '../data/coverage/{gsm}.transcripts.log'
    conda: "R"
    shell:
        """
        Rscript calculate_transcripts_coverage.R {output} {input.trim_distances} {params.study} {input.alignments} 2> {log}
        """


rule get_transcript_counts:
    """Sum the coverage within each transcript"""
    input:
        '../data/coverage/{gsm}.transcripts.bedGraph.gz'
    output:
        '../data/gene_counts/{gsm}.transcripts.tsv.gz'
    log:
        '../data/gene_counts/{gsm}.transcripts.log'
    conda: "R"
    shell:
        """
        Rscript calculate_transcript_counts.R {input} {output} 2> {log}
        """


rule get_codon_ribosome_density:
    """Get the average ribosome density around each codon"""
    params:
        transcript_annotations = lambda wildcards: merged_annotations.loc[merged_annotations['gsm'] == wildcards.gsm, 'transcript_annotations'].tolist()[0],
        motif_annotations = lambda wildcards: merged_annotations.loc[merged_annotations['gsm'] == wildcards.gsm, 'codon_annotations'].tolist()[0]
    input:
        '../data/coverage/{gsm}.transcripts.bedGraph.gz'
    output:
        motif_density_file = '../data/codon_density/{gsm}.tsv.gz',
        average_motif_density_file = '../data/codon_density/{gsm}.average.tsv.gz',
    log:
        '../data/codon_density/{gsm}.log'
    shell:
        """
        python find_ribosome_density_around_motifs.py \
          --read_count_cutoff 100 \
          --read_density_cutoff 0.33 \
          --cds_length_cutoff 100 \
          --l_overhang 50 \
          --r_overhang 50 \
          --transcript_annotation_file {params.transcript_annotations} \
          --motif_annotation_file {params.motif_annotations} \
          --ribosome_density_file {input} \
          --motif_density_file {output.motif_density_file} \
          --average_motif_density_file {output.average_motif_density_file} \
            &> {log}
        """


rule get_vk_motif_density:
    """Get the average ribosome density around each destabilizing amphipathic motif"""
    params:
        transcript_annotations = lambda wildcards: merged_annotations.loc[merged_annotations['gsm'] == wildcards.gsm, 'transcript_annotations'].tolist()[0],
        motif_annotations = lambda wildcards: merged_annotations.loc[merged_annotations['gsm'] == wildcards.gsm, 'vk_type_motif_annotations'].tolist()[0]
    input:
        '../data/coverage/{gsm}.transcripts.bedGraph.gz'
    output:
        motif_density_file = '../data/vk_type_density/{gsm}.tsv.gz'
    log:
        '../data/vk_type_density/{gsm}.log'
    shell:
        """
        python find_ribosome_density_around_motifs.py \
          --read_count_cutoff 100 \
          --read_density_cutoff 0.33 \
          --cds_length_cutoff 100 \
          --l_overhang 50 \
          --r_overhang 50 \
          --transcript_annotation_file {params.transcript_annotations} \
          --motif_annotation_file {params.motif_annotations} \
          --ribosome_density_file {input} \
          --motif_density_file {output.motif_density_file} \
            &> {log}
        """


rule plot_vk_motif_density:
    """Plot the average ribosome density around each destabilizing amphipathic motif"""
    input:
        motif_density_file = [f'../data/vk_type_density/{gsm}.tsv.gz' for gsm in set(merged_annotations.loc[:, "gsm"])]
    output:
        "plot_vk_type_density_han2020.nbconvert.ipynb"
    params:
        notebook = "plot_vk_type_density_han2020.ipynb"
    log:
        '../data/vk_type_density/plot_vk_type_density.log'
    conda: "R"
    shell:
        """
        export JUPYTER_DATA_DIR=$(pwd)
        export JUPYTER_CONFIG_DIR=$(pwd)
        jupyter nbconvert --to notebook --execute \
            --ExecutePreprocessor.kernel_name=ir {params.notebook} &> {log}
        """