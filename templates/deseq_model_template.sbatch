#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=20
#SBATCH -J libraryprotocol_librarydate_genotype_model
#SBATCH -o /scratch/mblab/chasem/rnaseq_pipeline/sbatch_log/libraryprotocol_librarydate_genotype_model_with_results.out

ml rnaseq_pipeline

raw_counts=/scratch/mblab/chasem/misc_scripts/deseq_model/data/exon_filtered_counts_20200812.csv
metadata=/scratch/mblab/chasem/misc_scripts/deseq_model/data/filtered_metadata_20200812.csv

output_dir=/scratch/mblab/chasem/misc_scripts/deseq_model/results


/scratch/mblab/chasem/misc_scripts/deseq_model/deseq_model.R -r ${raw_counts} -m ${metadata} -d ~LIBRARYPROTOCOL+LIBRARYDATE+GENOTYPE -g -o ${output_dir} -n libraryprotocol_librarydate_genotype_model_with_results
