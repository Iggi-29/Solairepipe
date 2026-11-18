snakemake -j20 -p --skip-script-cleanup --resources mem_mb=160000 --use-singularity --singularity-args "--cleanenv --bind /imppc/labs/eclab/" --use-conda --conda-frontend conda
