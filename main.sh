#!/bin/bash
#SBATCH --account=PAS1568
#SBATCH --time=1:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm-nf-rnaseq-%j.out
set -euo pipefail

# Load the Nextflow conda encironment
module load miniconda3/24.1.2-py310
conda activate nextflow

#Run the workflow
WORKFLOW=main.nf #eventually this should point to the github repository
nextflow run $WORKFLOW \
    -ansi-log false \
    -resume \
    "$@"
