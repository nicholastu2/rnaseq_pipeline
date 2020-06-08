#!/usr/bin/env nextflow

params.fastq_file_list = '/scratch/mblab/chasem/rnaseq_pipeline/job_scripts/nextflow_fastqfile_list_20200608_164325.csv' // this will need to be inputted by user

// split columns/rows of fastq_file_list for processing
Channel
    .fromPath(params.fastq_file_list)
    .splitCsv(header:true)
    .map{ row-> tuple(row.containing_directory, row.fastq_filename, file(row.organism), file(row.strandedness), file(row.output_dir)) }
    .groupTuple()
    .set { samples_channel }

process fastqc {
    input:
        set run_directory, fastq_filepath, organism, strandedness from samples_channel


    script:
    """
    echo $run_directory $fastq_filepath $organism $strandedness >> /scratch/mblab/chasem/nextflow_output_tester.txt
    """
}
