#!/bin/bash

#SBATCH --array=1-11%11
#SBATCH --mem=3G
#SBATCH -J raw_count
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/raw_count_20200924.out

ml rnaseq_pipeline

read count_dir < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/rnaseq_pipeline/job_scripts/env_pert_count_lookup_20200924.txt )

query_sheet=/scratch/mblab/chasem/rnaseq_pipeline/query/filtered_na_fastqFile_env_pert_20200923.csv

raw_count.py -c ${count_dir} -qs ${query_sheet} --interactive
