#!/bin/bash
#SBATCH -J Julia_genome_wget                 # A single job name for the array
#SBATCH -n 4                       # Number of cores
#SBATCH -N 1                       # All cores on one machine
#SBATCH --mem 8000                # in MB
#SBATCH -t 0-01:00                  # Maximum execution time (D-HH:MM)
#SBATCH -o job_%A_%a.log        # Standard output
#SBATCH -e job_%A_%a.log        # Standard error
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --account=mcgetthm-ra-icase
#SBATCH --qos castles

set -e

# This is a script for downloading stuff from the internet
# In this case we're getting the genome etc, from ensembl

module purge; module load bluebear # this line is required
module load wget/1.20.1-GCCcore-8.3.0

wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz