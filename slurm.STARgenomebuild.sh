#!/bin/bash
#SBATCH -J STAR_genome_build         # A single job name for the array
#SBATCH -n 4                       # Number of cores
#SBATCH -N 1                       # All cores on one machine
#SBATCH --mem 35000                 # in MB
#SBATCH -t 0-05:00                  # Maximum execution time (D-HH:MM)
#SBATCH -o job_%A_%a.log        # Standard output
#SBATCH -e job_%A_%a.log        # Standard error
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --account=mcgetthm-ra-icase
#SBATCH --qos castles

BB_WORKDIR=$(mktemp -d /scratch/${USER}_${SLURM_JOBID}.XXXXXX)
export TMPDIR=${BB_WORKDIR}

set -e

module purge; module load bluebear

module load STAR/2.7.2b-GCC-8.3.0

STAR --runThreadN ${SLURM_NTASKS} \
	 --runMode genomeGenerate \
	 --genomeDir /rds/projects/m/mcgetthm-ra-icase/Human \
	 --genomeFastaFiles /rds/projects/m/mcgetthm-ra-icase/Human/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	 --sjdbGTFfile /rds/projects/m/mcgetthm-ra-icase/Human/Homo_sapiens.GRCh38.101.gtf \
	 --sjdbOverhang 49

test -d ${BB_WORKDIR} && /bin/cp -r ${BB_WORKDIR} ./
test -d ${BB_WORKDIR} && /bin/rm -rf ${BB_WORKDIR}