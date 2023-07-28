#!/bin/bash

#SBATCH --account=BISC020703
#SBATCH --job-name=extract_grooming
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:0:0
#SBATCH --mem=50G

source activate FortMyrmidon_13March2023
Rscript /user/home/bzniks/code/ScriptsR_FINAL/STEP3_detect_grooming_on_all_data_CLUSTER.R $SLURM_ARRAY_TASK_ID


