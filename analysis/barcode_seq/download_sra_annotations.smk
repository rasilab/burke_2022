container: "../../burke_2022_latest.sif"

rule download_sra_annotations:
  input:
    "download_sra_annotations.ipynb"
  output:
    "../../annotations/sra_annotations.tsv"
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input}
    """
