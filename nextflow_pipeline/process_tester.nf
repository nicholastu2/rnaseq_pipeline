#!/usr/bin/env nextflow

params.fastq_file_list = 'fastq_file_list.csv' // this will need to be inputted by user

// split columns/rows of fastq_file_list for processing
Channel
    .fromPath(params.fastq_file_list)
    .splitCsv(header:true)
    .map{ row-> tuple(row.containing_directory, row.fastq_filename, file(row.organism), file(row.strandedness), file(row.output_dir)) }
    .groupTuple()
    .set { samples_channel }

process fastqc {
    input:


    script:
    """
    echo "I\'m working!" >> /scratch/mblab/chasem/nextflow_output_tester.txt
    """
}
