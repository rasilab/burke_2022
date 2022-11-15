container: "../burke_2022_0.8.2.sif"

rule all:
  input:
    "nanoluc_rrl_transit/scripts/analyze_alpha_beta_rrl_assay.nbconvert.ipynb",
    "nanoluc_rrl_transit/scripts/analyze_vk8_rrl_assay.nbconvert.ipynb",
    "circular_dichroism/scripts/cd_analysis.nbconvert.ipynb",
    "barcode_seq/pool1_mrna/scripts/plot_8xdicodon_effects.nbconvert.ipynb",
    "barcode_seq/pool1_mrna/scripts/plot_3d_pi_bulkiness_mrna.nbconvert.ipynb",
    "barcode_seq/pool1_mrna/scripts/plot_3d_pi_bulkiness_mrna.nbconvert.ipynb"


rule analyze_alpha_beta_rrl_assay:
  input:
    "nanoluc_rrl_transit/scripts/analyze_alpha_beta_rrl_assay.ipynb"
  output:
    "nanoluc_rrl_transit/scripts/analyze_alpha_beta_rrl_assay.nbconvert.ipynb" 
  conda: "R"  
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input}
    """

rule analyze_vk8_rrl_assay:
  input:
    "nanoluc_rrl_transit/scripts/analyze_vk8_rrl_assay.ipynb"
  output:
    "nanoluc_rrl_transit/scripts/analyze_vk8_rrl_assay.nbconvert.ipynb"
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
    "circular_dichroism/scripts/cd_analysis.nbconvert.ipynb"
  conda: "R"  
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input}
    """

rule pool1_mrna_8xdicodon_effects:
  input:
    "barcode_seq/pool1_mrna/scripts/plot_8xdicodon_effects.ipynb"
  output:
    "barcode_seq/pool1_mrna/scripts/plot_8xdicodon_effects.nbconvert.ipynb"
  conda: "R"  
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input}
    """

rule pool1_mrna_dipeptide_structure_effects:
  input:
    "barcode_seq/pool1_mrna/scripts/plot_dipeptide_secondary_structure_effects.ipynb"
  output:
    "barcode_seq/pool1_mrna/scripts/plot_dipeptide_secondary_structure_effects.nbconvert.ipynb"
  conda: "R"  
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input}
    """

rule pool1_mrna_3d_pi_bulkiness_effects:
  input:
    "barcode_seq/pool1_mrna/scripts/plot_3d_pi_bulkiness_mrna.ipynb"
  output:
    "barcode_seq/pool1_mrna/scripts/plot_3d_pi_bulkiness_mrna.nbconvert.ipynb"
  conda: "R"  
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input}
    """