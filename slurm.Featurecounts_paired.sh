#!/bin/bash
#SBATCH -J FeatureCounts    # A single job name for the array
#SBATCH -n 4                       # Number of cores
#SBATCH -N 1                       # All cores on one machine
#SBATCH --mem 35000                 # in MB
#SBATCH -t 0-20:00                 # Maximum execution time (D-HH:MM)
#SBATCH -o job_%A_%a.log        # Standard output
#SBATCH -e job_%A_%a.log        # Standard error
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --account=mcgetthm-macrophage-fibroblast-rnaseq
#SBATCH --qos castles

set -e

module purge; module load bluebear

module load bear-apps/2020a
module load Subread/2.0.1-GCC-8.3.0


# Summarize paired-end reads and count fragments (instead of reads):

featureCounts \
-p \
-t exon \
-g gene_id \
-a /rds/projects/m/mcgetthm-ra-icase/Human/Homo_sapiens.GRCh38.101.gtf \
-o /rds/projects/m/mcgetthm-ra-icase/Fibroblast-endothelial/Analysis/Counts/Featurecounts.txt \
/rds/projects/m/mcgetthm-ra-icase/Fibroblast-endothelial/Analysis/Aligned/*.out.bam

# -a = name of annotation file
# -o = name of output file
# -t = type of GTF annotation, 'exon' by default
# -g attribute type in GTF annotation, 'gene_id' by default
# --primary = counts primary alignments only 