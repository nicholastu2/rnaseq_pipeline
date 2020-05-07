#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --array=1-24%24
#SBATCH -D ./
#SBATCH -o sbatch_log/mblab_rnaseq_%A_%a.out
#SBATCH -e sbatch_log/mblab_rnaseq_%A_%a.err
#SBATCH -J mblab_rnaseq

ml novoalign/3.07.00
ml samtools/1.6
ml htseq/0.9.1
read fastq_file < <( sed -n ${SLURM_ARRAY_TASK_ID}p job_scripts/run_0629_0618_fastq_list.txt ); set -e

mkdir -p reports/run_0629_0618
sample=${fastq_file##*/}; sample=${sample%.f*q.gz}; novoalign -c 8 -o SAM -d KN99/crNeoKN99.nix -f ${fastq_file} 2> reports/run_0629_0618/${sample}_novoalign.log | samtools view -bS > reports/run_0629_0618/${sample}_aligned_reads.bam
sample=${fastq_file##*/}; sample=${sample%.f*q.gz}; novosort --threads 8 reports/run_0629_0618/${sample}_aligned_reads.bam > reports/run_0629_0618/${sample}_sorted_aligned_reads.bam 2> reports/run_0629_0618/${sample}_novosort.log
sample=${fastq_file##*/}; sample=${sample%.f*q.gz}; htseq-count -f bam -s no -t CDS reports/run_0629_0618/${sample}_sorted_aligned_reads.bam KN99/crNeoKN99.gtf > reports/run_0629_0618/${sample}_read_count.tsv 2> reports/run_0629_0618/${sample}_htseq.log
