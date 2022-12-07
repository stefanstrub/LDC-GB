#!/usr/bin/env bash
#PBS -l walltime=71:00:00
#PBS -l select=1:ncpus=1:mem=4G
#PBS -N snk_main
module load python/3.7.2
module load fftw
module load singularity
source $HOME/venv/bin/activate
export SINGULARITY_BINDPATH="/work/SC/lisa/LDC/ancillary_data:/data"
cd /work/SC/lisa/lejeune/sangria
snakemake -j 100 --use-singularity --default-resources "mem_mb=4000" --cluster "qsub -l walltime=36:00:00,select=1:ncpus=1:mem={resources.mem_mb}mb" 
