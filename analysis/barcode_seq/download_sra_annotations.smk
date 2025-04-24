rule download_sra_annotations:
  input:
    "download_sra_annotations.ipynb"
  output:
    "../../annotations/sra_annotations.tsv"
  container: "docker://ghcr.io/rasilab/reutils_geoquery:1.0.0"
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input}
    """
