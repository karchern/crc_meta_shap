#!/bin/bash -e
#SBATCH -A zeller
#SBATCH --time=0-24:00:00
#SBATCH --qos=low
## Looks like we need more than 100 GB of ram for the analysis...
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
### The following number might need to be updated!
#SBATCH --array=1-100
#SBATCH --error=/scratch/karcher/shap/%a.err
#SBATCH --output=/scratch/karcher/shap/%a.out
#SBATCH --job-name=get_shap

cat iv.sh |  head -n $SLURM_ARRAY_TASK_ID | tail -1 | bash
