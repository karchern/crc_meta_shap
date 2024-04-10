#!/bin/bash -e
#SBATCH -A zeller
#SBATCH --time=0-04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --array=1-200
#SBATCH --error=/scratch/karcher/shap/%a.err
#SBATCH --output=/scratch/karcher/shap/%a.out
#SBATCH --job-name=get_shap

## RUN THIS FROM /g/scb/zeller/karcher/anaconda3/envs/r_test_new
cd /g/scb/zeller/karcher/SHAP/scripts
cat iv.sh |  head -n $SLURM_ARRAY_TASK_ID | tail -1 | bash
