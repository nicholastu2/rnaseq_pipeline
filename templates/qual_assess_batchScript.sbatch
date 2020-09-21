#!/bin/bash

#SBATCH --array=1-4%4
#SBATCH --mem=3G
#SBATCH -J qual_assess
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/env_pert_qa.out

ml rnaseq_pipeline

read run_dir < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/rnaseq_pipeline/experiments/qual_assess_env_lookup.txt )

query_sheet=/scratch/mblab/chasem/rnaseq_pipeline/query/environmental_perturbation_20200817.csv

echo "working on ${run_dir}"

quality_assess_1.py -ac ${run_dir} -qs ${query_sheet} --interactive
