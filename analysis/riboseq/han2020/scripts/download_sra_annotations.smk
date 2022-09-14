container:  '../../../../burke_2022_latest.sif'

rule all:
    """List of all files we want at the end"""
    input:
        '../annotations/sra_annotations.tsv',


rule get_sra_annotations:
    input:
        '../annotations/geo_accession_numbers.csv'
    output:
        '../annotations/sra_annotations.tsv'
    params:
        notebook = 'download_geo_sra_annotations_based_on_geo_accession_numbers.Rmd',
    conda: "R"
    shell:
        """
        Rscript -e "rmarkdown::render('{params.notebook}')"
        """