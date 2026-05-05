#!/bin/bash -e

#SBATCH --job-name=frag_sim
#SBATCH --mail-user=j.n.thiem@fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=700
#SBATCH --time=03:00:00
#SBATCH --qos=standard
#SBATCH --output=logs/frag_sim_%A_%a.out
#SBATCH --error=logs/frag_sim_%A_%a.err
#SBATCH --array=1-30  

module add R
module add GDAL

Rscript "extended_pipeline.R" ${SLURM_ARRAY_TASK_ID} ${SLURM_ARRAY_JOB_ID}
