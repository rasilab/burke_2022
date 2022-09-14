container:
    "../../../../burke_2022_latest.sif"


rule all:
    input:
        # contaminants fasta
        "../data/reference_sequences/hg38.rrna.fasta",
        # gencode annotations
        "../data/gencode/gencode.v32.annotation.gff3.gz",
        # gencode unique CCDS transcripts
        "../data/gencode/gencode.v32.canonical_ccds_tx.fa",
        # bowtie2 contaminant index
        "../data/bowtie2/hg38.rrna.1.bt2",
        # bowtie2 transcript database
        "../data/gencode/hg38.gencode.v32.canonical_ccds_tx.1.bt2",


rule download_contaminants:
    input:
        notebook = "download_rrna_contaminant_sequences.ipynb"
    output:
        "../data/reference_sequences/hg38.rrna.fasta"
    shell:
        """
        export JUPYTER_DATA_DIR=$(pwd)
        export JUPYTER_CONFIG_DIR=$(pwd)
        jupyter nbconvert --to notebook --execute {input}
        """


rule download_gencode_annotations:
    output:
        "../data/gencode/gencode.v32.annotation.gff3.gz"
    shell:
        """
        wget \
        ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gff3.gz \
        -O ../data/gencode/gencode.v32.annotation.gff3.gz
        """


rule parse_unique_ccds_transcripts:
    input:
        annotation = "../data/gencode/gencode.v32.annotation.gff3.gz"
    output:
        fasta = "../data/gencode/gencode.v32.canonical_ccds_tx.fa",
        annotation = "../data/gencode/gencode.v32.canonical_ccds.parameters.tsv.gz"
    conda: "R"
    shell:
        """
        Rscript make_unique_ccds_transcripts.R
        """


rule make_contaminant_bowtie_database:
    input: "../data/reference_sequences/hg38.rrna.fasta"
    output:  "../data/bowtie2/hg38.rrna.1.bt2"
    params:  
        bowtie2_prefix = "../data/bowtie2/hg38.rrna"
    conda: "bowtie2_samtools"
    threads: 36
    shell:
        """
        bowtie2-build --threads {threads} {input} {params.bowtie2_prefix}
        """


rule make_ccds_transcript_database:
    input: "../data/gencode/gencode.v32.canonical_ccds_tx.fa"
    output: "../data/gencode/hg38.gencode.v32.canonical_ccds_tx.1.bt2"
    params:
        bowtie2_prefix = "../data/gencode/hg38.gencode.v32.canonical_ccds_tx"
    conda: "bowtie2_samtools"
    threads: 36
    shell:
        """
        bowtie2-build --threads {threads} {input} {params.bowtie2_prefix}
        """
