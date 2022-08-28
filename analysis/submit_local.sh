snakemake \
    --use-conda \
    --cores=all \
    --use-singularity \
    --singularity-args "--bind /fh" \
    -p $@
