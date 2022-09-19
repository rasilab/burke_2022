container:
    "../../../../burke_2022_latest.sif"

rule all:
    input:
        ['../data/motif_counts/ccds_codon_counts_frame' + str(n) + '.tsv' for n in range(3)],
        '../data/motif_counts/ccds_codon_locs.tsv.gz',
        '../data/motif_counts/destabilizing_vk_motifs.tsv.gz',
        '../data/motif_counts/kk_motifs.tsv.gz',


rule get_cds_codon_freq:
    """Get the frequency of codons in CDS"""
    input:
        notebook = "calculate_ccds_codon_freq.ipynb"
    output:
        codon_motif_file = ['../data/motif_counts/ccds_codon_counts_frame' + str(n) + '.tsv' for n in range(3)],
    params:
        transcript_seq = '../data/gencode/gencode.v32.canonical_ccds.fa',
    log:
        '../data/motif_counts/codon_counts.log'
    conda:
        'base'
    shell:
        """
        export JUPYTER_DATA_DIR=$(pwd)
        export JUPYTER_CONFIG_DIR=$(pwd)
        jupyter nbconvert --to notebook --execute {input.notebook}
        """


rule get_cds_codon_loc:
    """Get the location of codons in CDS"""
    input:
        notebook = "find_codon_locations_in_ccds.ipynb"
    output:
        "../data/motif_counts/ccds_codon_locs.tsv.gz"

    log:
        '../data/motif_counts/codon_locs.log'
    conda:
        'base'
    shell:
        """
        export JUPYTER_DATA_DIR=$(pwd)
        export JUPYTER_CONFIG_DIR=$(pwd)
        jupyter nbconvert --to notebook --execute {input.notebook}
        """


rule get_destabilizing_motif_counts:
    """Get the frequency of destabilizing in CDS"""
    input:
        notebook = "find_destabilizing_vk_motifs.ipynb"
    output:
        motif_file = '../data/motif_counts/destabilizing_vk_motifs.tsv.gz',
    params:
        transcript_seq = '../data/gencode/gencode.v32.canonical_ccds.fa',
    log:
        '../data/motif_counts/destabilizing_vk_motifs_counts.log'
    conda:
        'base'
    shell:
        """
        export JUPYTER_DATA_DIR=$(pwd)
        export JUPYTER_CONFIG_DIR=$(pwd)
        jupyter nbconvert --to notebook --execute {input.notebook}
        """


rule get_kk_motif_counts:
    """Get the frequency of destabilizing in CDS"""
    input:
        notebook = "find_kk_motifs.ipynb"
    output:
        motif_file = '../data/motif_counts/kk_motifs.tsv.gz',
    params:
        transcript_seq = '../data/gencode/gencode.v32.canonical_ccds.fa',
    log:
        '../data/motif_counts/kk_motifs_counts.log'
    conda:
        'base'
    shell:
        """
        export JUPYTER_DATA_DIR=$(pwd)
        export JUPYTER_CONFIG_DIR=$(pwd)
        jupyter nbconvert --to notebook --execute {input.notebook}
        """
        