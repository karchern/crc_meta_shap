#!/bin/bash -e
#SBATCH -A zeller
#SBATCH --time=0-48:00:00
## Looks like we need more than 100 GB of ram for the analysis...
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
### The following number might need to be updated!
#SBATCH --array=1-150
#SBATCH --error=/scratch/karcher/shap/%a.err
#SBATCH --output=/scratch/karcher/shap/%a.out
#SBATCH --job-name=get_shap

## RUN THIS FROM conda activate /g/scb2/zeller/karcher/mambaforge/envs/mlr3_0.20.2
cd /g/scb/zeller/karcher/SHAP/scripts
cat iv.sh |  head -n $SLURM_ARRAY_TASK_ID | tail -1 | bash
