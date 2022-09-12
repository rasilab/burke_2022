container:
    "../../../../burke_2022_latest.sif"

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