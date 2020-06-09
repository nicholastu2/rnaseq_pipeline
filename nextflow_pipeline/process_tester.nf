#!/usr/bin/env nextflow

params.fastq_file_list = '/scratch/mblab/chasem/rnaseq_pipeline/job_scripts/nextflow_fastqfile_list_20200608_164325.csv' // this will need to be inputted by user
params.scratch_sequence = '/scratch/mblab/mblab.shared/scratch_sequence'

// split columns/rows of fastq_file_list for processing
Channel
    .fromPath(params.fastq_file_list)
    .splitCsv(header:true)
    .map{ row-> tuple(row.runDirectory, file(row.fastqFileName), row.organism, row.strandedness) }
    .set { samples_channel }

scratch_sequence = file(params.scratch_sequence)

// this works
process make_reports_directory {
    cache = 'false'
    executor = 'local'

    input:
        set val(run_directory), file(fastq_filepath), val(organism), val(strandedness) from samples_channel
    output:
        tuple run_directory, organism, strandedness, fastq_filepath into samples_channel_ch

    """
    mkdir -p ${scratch_sequence}/${run_directory}
    echo $organism >> /scratch/mblab/chasem/nextflow_output_tester.txt
    """
}

process novoalign {
    stageInMode = 'copy'
    stageOutMode = 'rsync'
    beforeScript = 'module load novoalign; module load samtools'
    cache = 'false'
    executor = 'slurm'
    publishDir "/scratch/mblab/chasem/rnaseq_pipeline/reports", mode: 'copy'

    input:
        set run_directory, organism, strandedness, fastq_filepath from samples_channel_ch

    script:
    fastq_simple_name = fastq_filepath.getSimpleName()
    """
    echo "novoalign -r All -c 8 -o SAM -d genome_files -r all -f ${fastq_simple_name} 2> ${fastq_simple_name}_novoalign.log \
    | samtools view -bS > ${fastq_simple_name}_aligned_reads.bam" >> /scratch/mblab/chasem/nextflow_output_tester.txt

    """
}