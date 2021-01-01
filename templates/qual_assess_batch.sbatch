#!/bin/bash

#SBATCH --array=1-56%30
#SBATCH --mem=10G
#SBATCH -J qual_assess
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/qual_assess_1.out

ml rnaseq_pipeline

read run_dir < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/rnaseq_pipeline/job_scripts/induction_set_qual_assess_lookup.txt )

query=/scratch/mblab/chasem/rnaseq_pipeline/query/induction_set_query.csv

quality_assess_1.py -ac ${run_dir} -qs ${query} --interactive
