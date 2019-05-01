#!/bin/bash -l
#SBATCH --time=0-00:05:00 --mem-per-cpu=500
#SBATCH -p debug
#SBATCH -o job-%a.out
#SBATCH --array=0-9
module load matlab
matlab -nojvm -r "simulate_GMS_Ju($SLURM_ARRAY_TASK_ID,16,100,1,0); quit"