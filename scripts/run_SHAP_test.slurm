#!/bin/bash -e
#SBATCH -A zeller
#SBATCH --time=0-72:00:00
## Looks like we need more than 100 GB of ram for the analysis...
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
### The following number might need to be updated!
#SBATCH --array=1-50
#SBATCH --error=/scratch/karcher/shap/%a.err
#SBATCH --output=/scratch/karcher/shap/%a.out
#SBATCH --job-name=get_shap

## RUN THIS FROM /g/scb/zeller/karcher/anaconda3/envs/r_test_new


cd /g/scb/zeller/karcher/SHAP/scripts
#cat iv.sh |  head -n $SLURM_ARRAY_TASK_ID | tail -1 | bash
cat iv_training.sh |  head -n $SLURM_ARRAY_TASK_ID | tail -1 | bash
