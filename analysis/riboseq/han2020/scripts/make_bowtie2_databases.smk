container:
    "../../../../burke_2022_latest.sif"

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