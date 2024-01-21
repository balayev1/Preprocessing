#!/bin/bash
#SBATCH -p pi_kaminski,day --time=8:00:00 --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=32GB  --job-name=demultiplexing -o /vast/palmer/scratch/kaminski/ab3478/Demultiplex/demuplexing.out -e /vast/palmer/scratch/kaminski/ab3478/Demultiplex/demuplexing.err

module load bcl2fastq2/2-20-0-foss-2018b

#######################################
################################
### Demultiplexing snRNA-seq data
/gpfs/gibbs/pi/kaminski/public/softwares/cellranger-4.0.0/cellranger-4.0.0/bin/cellranger mkfastq --run /gpfs/ysm/pi/haifan_lin/sequencers/runs/220922_K00175_0283_BHT5TFBBXY \
    --lanes=1,2 \
    --simple-csv=/vast/palmer/scratch/kaminski/ab3478/Demultiplex/demultiplex_1052022.csv \
    --localcores=8 \
    --mempercore=25 \
    --force-single-index \
    --output-dir=/vast/palmer/scratch/kaminski/ab3478/snRNA_forATAC \
    --use-bases-mask=Y28n*,I10n*,N10,Y90n* 
