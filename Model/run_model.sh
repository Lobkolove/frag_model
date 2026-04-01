#!/bin/bash -e

#SBATCH --job-name=frag_red_hab
#SBATCH --mail-user=j.n.thiem@fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=700
#SBATCH --time=00:05:00
#SBATCH --qos=standard

module add R
module add GDAL

Rscript "sim_pipeline.R"
