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
        set run_directory, file(fastq_filepath), organism, strandedness from samples_channel
    output:
        tuple(run_directory, organism, strandedness, fastq_filepath) into samples_channel_ch

    """
    mkdir -p ${scratch_sequence}/${run_directory}
    """
}

process novoalign {
    stageInMode = 'copy'
    stageOutMode = 'rsync'
    beforeScript = 'module load novoalign'
    cache = 'false'
    executor = 'sge'
    publishDir ${scratch_sequence}/${run_directory}, mode: 'copy'

    input:
        set run_directory, organism, strandedness, fastq_filepath from samples_channel_ch
    output:
        tuple(novoalign_output, novoalign_log) into novoalign_ch

    script:
    """
    echo "novoalign -r All -c 8 -o SAM -d param.KN99_genome_index -f ${fastq_filepath} 2> ${organism}_novoalign.log)" >> /scratch/mblab/chasem/nextflow_output_tester.txt

    """
}
novoalign_ch.subscribe{ println + it.novoalign_log }
