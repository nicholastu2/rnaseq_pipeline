#!/usr/bin/env nextflow

params.fastq_file_list = '/scratch/mblab/chasem/rnaseq_pipeline/job_scripts/nextflow_fastqfile_list_20200608_164325.csv' // this will need to be inputted by user
params.scratch_sequence = '/scratch/mblab/mblab.shared/scratch_sequence'

// split columns/rows of fastq_file_list for processing
Channel
    .fromPath(params.fastq_file_list)
    .splitCsv(header:true)
    .map{ row-> tuple(row.runDirectory, row.fastqFileName, row.organism, row.strandedness) }
    .groupTuple()
    .set { samples_channel }


process lts_to_scratch {
    cache = 'false'
    executor = 'local'
    beforeScript = 'ml rnaseq_pipeline'

    input:
        set run_directory, fastq_filepath, organism, strandedness from samples_channel

    // extract parent directory from source, make directory of same name to deposit file in destination
    script:
    """
    for fastqfile in ${fastq_filepath}
    do
        echo "rsync -aHv --sample ${fastqfile} -d  param.scratch_sequence/${run_directory}" >  /scratch/mblab/chasem/nextflow_output_tester.txt //somehow this needs to go to a run_#### file in scratch sequence b/c fastqfilenames may not be unique
    done
    """
}