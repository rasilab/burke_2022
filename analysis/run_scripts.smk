container: "burke_2022_latest.sif"

rule analyze_alpha_beta_rrl_assay:
  input:
    "nanoluc_rrl_transit/scripts/analyze_alpha_beta_rrl_assay.ipynb"
  conda: "R"  
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input}
    """

rule cd_analysis:
  input:
    "circular_dichroism/scripts/cd_analysis.ipynb"
  output:
    "circular_dichroism/figures/sesca_ss_content.pdf"
    "circular_dichroism/figures/cd_spectra.pdf"
  conda: "R"  
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input}
    """
