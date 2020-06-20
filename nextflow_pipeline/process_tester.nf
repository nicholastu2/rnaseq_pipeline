#!/usr/bin/env nextflow

// split columns/rows of fastq_file_list for processing
Channel
    .fromPath(params.fastq_file_list)
    .splitCsv(header:true)
    .map{row-> tuple(row.runDirectory, file(row.fastqFileName), row.organism, row.strandedness) }
    .set { fastq_filelist }


scratch_sequence = file(params.scratch_sequence)

// move files out of /lts into scratch work directory (make available for slurm)
process toScratch {
    // errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    // disk
    executor "local"
    stageInMode "copy"
    stageOutMode "move"

    input:
        tuple val(run_directory), file(fastq_filepath), val(organism), val(strandedness) from fastq_filelist
    output:
        tuple val(run_directory), file(fastq_filepath), val(organism), val(strandedness) into fastq_filelist_ch
        tuple val(run_directory), file(fastq_filepath) into fastqc_ch

    script:
        """
        """
}

fastq_ch_tuples = fastqc_ch.collect().groupTuple(by: run_directory)

process fastqc {
    echo true
    executor "slurm"
    memory "12G"
    cpus 8
    beforeScript "ml fastqc"
    publishDir "$params.align_count_results/$run_directory/fastqc", mode:"copy", overwite: true

    input:
        set run_number, reads from fastq_ch_tuples

    output:
        file "*_fastqc.{zip,html}" into fastqc_results

    script:
      """
      fastqc --quiet --threads 8 ${reads}
      """
}

// process files in work directory with slurm
process novoalign {
    scratch true
    executor "slurm"
    cpus 8
    memory "12G"
    beforeScript "ml novoalign"
    stageInMode "copy"
    stageOutMode "move"
    afterScript "rm ${fastq_file}" // figure out how to actually delete these


    input:
        tuple val(run_directory), file(fastq_file), val(organism), val(strandedness) from fastq_filelist_ch
    output:
        tuple val(run_directory), val(fastq_simple_name), val(organism), val(strandedness), file("${fastq_simple_name}_aligned_reads.sam") into sam_align_ch
        tuple val(run_directory), file("${fastq_simple_name}_novoalign.log") into novoalign_log_ch

    script:
        fastq_simple_name = fastq_file.getSimpleName()
        if (organism == 'S288C_R64')
            """
            novoalign -r All -c 8 -o SAM -d ${params.S288C_R64_novoalign_index} -f ${fastq_file} 1> ${fastq_simple_name}_aligned_reads.sam 2> ${fastq_simple_name}_novoalign.log
            """
        else if (organism == 'KN99')
            """
            novoalign -r All -c 8 -o SAM -d ${params.KN99_novoalign_index} -f ${fastq_file} 1> ${fastq_simple_name}_aligned_reads.sam 2> ${fastq_simple_name}_novoalign.log
            """
}

process novosort {
  scratch true
  executor "slurm"
  memory "12G"
  cpus 8
  stageInMode "copy"
  stageOutMode "move"
  beforeScript "ml novoalign"
  publishDir "$params.align_count_results/$run_directory", mode:"move", overwite: true, pattern: "*_novosort.log"

    input:
      tuple val(run_directory), val(fastq_simple_name), val(organism), val(strandedness), file(alignment_sam) from sam_align_ch
    output:
      tuple val(run_directory), val(fastq_simple_name), val(organism), val(strandedness), file("${fastq_simple_name}_sorted_aligned_reads.bam") into sorted_bam_align_ch
      file "${fastq_simple_name}_novosort.log" into novosort_ch

    // check what happens with the --index option
    // see http://www.novocraft.com/documentation/novosort-2/ for novosort options

    script:
      """
      samtools view -bS ${alignment_sam} | novosort --threads 8 --markDuplicates --index --output ${fastq_simple_name}_sorted_aligned_reads.bam 2> ${fastq_simple_name}_novosort.log
      """
}
