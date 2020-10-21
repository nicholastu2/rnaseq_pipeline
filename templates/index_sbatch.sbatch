#!/bin/bash

#SBATCH --array=1-1974%30
#SBATCH --mem=3G
#SBATCH -J index_align_count
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/index_align_count

ml rnaseq_pipeline

read bam_file < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/rnaseq_pipeline/bam_lookup.txt )

samtools index ${bam_file}
