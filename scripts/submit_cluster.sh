snakemake \
    --jobs 999 \
    --cluster-config cluster.yaml \
    --cluster "sbatch -n {cluster.n}  -t {cluster.time}" \
    --use-conda \
    --use-singularity \
    --singularity-args "--bind /fh:/fh" \
    --cores=all \
    -p
