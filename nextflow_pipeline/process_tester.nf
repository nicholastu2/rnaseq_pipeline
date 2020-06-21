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
        //tuple val(run_directory), file(fastq_filepath) into fastqc_ch

    script:
        """
        """
}

fastq_filelist_ch.into { fastq_filelist_ch; fastqc }

process fastqc {

    scratch true
    executor "slurm"
    memory "12G"
    cpus 8
    beforeScript "ml fastqc"
    publishDir "$params.align_count_results/$run_directory/fastqc", mode:"copy", overwite: true, pattern: "*_fastqc.{zip,html}"

    input:
        tuple val(run_directory), values from fastqc_ch.collect().groupBy()

    output:
        file "*_fastqc.{zip,html}"

    script:
      flat_read_list = $values[0].flatten()
      """
      fastqc --quiet --threads 8 ${flat_read_list}
      """
}

// process files in work directory with slurm
process novoalign {

    scratch true
    executor "slurm"
    cpus 8
    memory "12G"
    beforeScript "ml novoalign samtools"
    stageOutMode "move"
    // afterScript "rm ${fastq_file}" // figure out how to actually delete these
    publishDir "$params.align_count_results/$run_directory/logs", mode:"copy", overwite: true, pattern: "*_novo*.log"
    publishDir "$params.align_count_results/$run_directory/align", mode:"copy", overwite: true, pattern: "*_reads.bam.*"


    input:
        tuple val(run_directory), file(fastq_file), val(organism), val(strandedness) from fastq_filelist_ch
    output:
        tuple val(run_directory), val(fastq_simple_name), val(organism), val(strandedness), file("${fastq_simple_name}_sorted_aligned_reads.bam") into bam_align_ch
        tuple val(run_directory), file("${fastq_simple_name}_novoalign.log") into novoalign_log_ch
        file("${fastq_simple_name}_novosort.log") into novosort_log_ch

    script:
        fastq_simple_name = fastq_file.getSimpleName()
        if (organism == 'S288C_R64')
            """
            novoalign -r All -c 8 -o SAM -d ${params.S288C_R64_novoalign_index} -f ${fastq_file} 2> ${fastq_simple_name}_novoalign.log | samtools view -bS | novosort - --threads 8 --markDuplicates -o ${fastq_simple_name}_sorted_aligned_reads.bam 2> ${fastq_simple_name}_novosort.log
            """
        else if (organism == 'KN99')
            """
            novoalign -r All -c 8 -o SAM -d ${params.KN99_novoalign_index} -f ${fastq_file} 2> ${fastq_simple_name}_novoalign.log | samtools view -bS | novosort - --threads 8 --markDuplicates -o ${fastq_simple_name}_sorted_aligned_reads.bam 2> ${fastq_simple_name}_novosort.log
            """
}

process htseqCount {

  scratch true
  executor "slurm"
  memory "12G"
  beforeScript "ml htseq"
  publishDir "$params.align_count_results/$run_directory/count", mode:"copy", overwite: true, pattern: "*_read_count.tsv"
  publishDir "$params.align_count_results/$run_directory/log", mode:"copy", overwite: true, pattern: "*_htseq.log"

    input:
      tuple val(run_directory), val(fastq_simple_name), val(organism), val(strandedness), file(sorted_alignment_bam) from bam_align_ch

    output:
      tuple val(run_directory), file(sorted_alignment_bam), file("${fastq_simple_name}_read_count.tsv") into align_count_output_ch
      file("${fastq_simple_name}_htseq.log") into htseq_log_ch
      file("${fastq_simple_name}_htseq_annote.sam") into htseq_sam_ch

    script:
    // note: this is set for gff3 format. eg: (tab delimited)
    // CP022321.1	EuPathDB	exon	1767728	1768158	.	+	.	ID=exon_CKF44_00681-E1;Parent=CKF44_00681-t42_2,CKF44_00681-t42_1
    // hence -t exon -i ID
      if (organism == 'KN99')
        """
        htseq-count -f bam -o ${fastq_simple_name}_htseq_annote.sam -s ${strandedness} -t exon ${sorted_alignment_bam} ${params.KN99_annotation_file} 1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log
        """
      else if (organism == 'S288C_R64')
        """
        htseq-count -f bam -o ${fastq_simple_name}_htseq_annote.sam -s ${strandedness} -t gene -i ID ${sorted_alignment_bam} ${params.S288C_R64_annotation_file} 1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log
        """
      else if (organism == 'H99')
        """
        htseq-count -f bam -s ${strandedness} -t exon -i gene_id ${sorted_alignment_bam} ${params.H99_annotation_file} 1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log
        """
      else
        error "Invalid Organism"

}
