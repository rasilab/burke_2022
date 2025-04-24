snakemake \
    --cores=all \
    --use-singularity \
    --singularity-args "--bind /fh" \
    -p $@
