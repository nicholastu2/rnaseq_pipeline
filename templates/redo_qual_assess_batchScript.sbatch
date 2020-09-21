#!/bin/bash

#SBATCH --array=1-9%9
#SBATCH --mem=2G
#SBATCH -J qual_assess
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/redo_failed_qual_assess.out

ml rnaseq_pipeline

read run_dir < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/misc_scripts/failed_qual_assess_lookup.txt )

quality_assess_1.py -ac ${run_dir} -o ${run_dir} -qs /scratch/mblab/chasem/rnaseq_pipeline/query/nextflow_full_test.csv --interactive
