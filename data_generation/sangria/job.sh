#!/usr/bin/env bash
#$ -P P_lisaf
#$ -l h_rt=172:00:00
#$ -l s_rss=4G
#$ -l sps=1
#$ -N snk_main
#$ -j y
#$ -q longlasting
##$ -l demon=1
source /sps/hep/lisaf/mlejeune/venv/bin/activate
cd /sps/hep/lisaf/mlejeune/sangria/gb-pipeline
ccenv singularity 3.5.2
export SINGULARITY_BINDPATH="/sps/hep/lisaf/mlejeune/data:/data"
#snakemake -j 1 --use-singularity --cluster "ssh mlejeune@cca004 /opt/sge/bin/lx-amd64/qsub -l h_rt=00:10:00,s_rss=1G -P P_lisaf"
snakemake -j 1000 --use-singularity --cluster "ssh mlejeune@cca004 /opt/sge/bin/lx-amd64/qsub -l h_rt=48:00:00,s_rss=4G -P P_lisaf -l sps=1" --batch batch_merge=2/2
