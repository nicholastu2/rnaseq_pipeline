#!/bin/bash

#SBATCH --array=1-56%30
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=20
#SBATCH -J deseq
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/deseq.out

ml rnaseq_pipeline

read dds < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/rnaseq_pipeline/job_scripts/dds_lookup.txt )

output_dir=/scratch/mblab/chasem/rnaseq_pipeline/experiments/ninety_min_induct

deseq_model.R -d ${dds} -o ${output_dir} -n full_model
