#!/bin/bash
#SBATCH -J Multiqc      # A single job name for the array
#SBATCH -n 4                       # Number of cores
#SBATCH -N 1                       # All cores on one machine
#SBATCH --mem 8000                # in MB
#SBATCH -t 1-00:00                  # Maximum execution time (D-HH:MM)
#SBATCH -o job_%A_%a.log        # Standard output
#SBATCH -e job_%A_%a.log        # Standard error
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --account=mcgetthm-macrophage-fibroblast-rnaseq
#SBATCH --qos castles

set -e
module purge; module load bluebear # this line is required

module load MultiQC/1.9-foss-2019b-Python-3.7.4

multiqc ./Fasta_files
multiqc ./Aligned
multiqc ./Counts