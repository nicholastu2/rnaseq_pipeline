#!/bin/bash

#SBATCH --array=1-61%29
#SBATCH --mem=20G
#SBATCH -J raw_count
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/raw_count.out

ml rnaseq_pipeline

read count_dir < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/misc_scripts/count_directory_lookup.txt )

query_sheet=/scratch/mblab/chasem/rnaseq_pipeline/query/combined_df_20200817.csv

raw_count.py -c ${count_dir} -qs ${query_sheet} --interactive
