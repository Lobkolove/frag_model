#!/bin/bash -e

#SBATCH --job-name=frag_seq_test
#SBATCH --mail-user=j.n.thiem@fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=02:00:00
#SBATCH --qos=standard
#SBATCH --output=logs/frag_seq_test_%j_%a.out
#SBATCH --array=1-15

module add R
module add GDAL

Rscript "sim_pipeline.R" ${SLURM_ARRAY_TASK_ID}
