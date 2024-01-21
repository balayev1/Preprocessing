#!/bin/bash
#SBATCH --time=48:00:00 --cpus-per-task=8 --mem=16000M  --job-name=featureCounts.PD1I -o /home/abalay/scratch/PD1I_datasets/Transcriptomics/featureCounts.PD1I.sh.o%J -e /home/abalay/scratch/PD1I_datasets/Transcriptomics/featureCounts.PD1I.sh.e%J

## Activate conda environment
module load anaconda3
source activate agshin_env

Rscript featureCounts_script.r

