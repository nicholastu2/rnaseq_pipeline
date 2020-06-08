#!/usr/bin/env nextflow

params.fastq_file_list = '/scratch/mblab/chasem/rnaseq_pipeline/job_scripts/nextflow_fastqfile_list_20200608_164325.csv' // this will need to be inputted by user
params.scratch_sequence = '/scratch/mblab/mblab.shared/scratch_sequence'

// split columns/rows of fastq_file_list for processing
Channel
    .fromPath(params.fastq_file_list)
    .splitCsv(header:true)
    .map{ row-> [row.runDirectory: row.fastqFileName, row.fastqFileName: row.organism, row.fastqFileName: row.strandedness] }
    .set { samples_channel }

scratch_sequence = file(params.scratch_sequence)

process make_scratch_directory {
    cache = 'false'
    executor = 'local'

    input:
        set row_array from samples_channel

    script:
        """
        echo $row_array

        """
}

