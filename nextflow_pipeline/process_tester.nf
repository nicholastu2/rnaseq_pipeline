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

    script:
        """
        """
}

// process files in work directory with slurm
process alignCount {

    scratch true
    executor "slurm"
    cpus 8
    memory "12G"
    beforeScript "ml novoalign samtools htseq"
    stageOutMode "move"
    // afterScript "rm ${fastq_file}" // figure out how to actually delete these
    publishDir "$params.align_count_results/$run_directory/logs", mode:"copy", overwite: true, pattern: "*.log"
    publishDir "$params.align_count_results/$run_directory/count", mode:"copy", overwite: true, pattern: "*_read_count.tsv"
    publishDir "$params.align_count_results/$run_directory/align", mode:"copy", overwite: true, pattern: "*sorted_aligned_reads_reads_with_annote.bam.*"


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
            novoalign -r All -c 8 -o SAM -d ${params.S288C_R64_novoalign_index} -f ${fastq_file} 2> ${fastq_simple_name}_novoalign.log | \
            samtools view -bS |
            novosort - --threads 8 --markDuplicates -o ${fastq_simple_name}_sorted_aligned_reads.bam 2> ${fastq_simple_name}_novosort.log

            htseq-count -f bam -o ${fastq_simple_name}_htseq_annote.sam -s ${strandedness} -t exon ${fastq_simple_name}_sorted_aligned_reads.bam ${params.KN99_annotation_file} 1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log

            sed "s/\t//" ${fastq_simple_name}_htseq_annote.sam > ${fastq_simple_name}_htseq_annote.sam

            samtools view ${fastq_simple_name}_sorted_aligned_reads.bam | paste - b > ${fastq_simple_name}_sorted_aligned_reads_with_annote.bam
            """
        else if (organism == 'KN99')
            """
            novoalign -r All -c 8 -o SAM -d ${params.KN99_novoalign_index} -f ${fastq_file} 2> ${fastq_simple_name}_novoalign.log | \
            samtools view -bS | \
            novosort - --threads 8 --markDuplicates -o ${fastq_simple_name}_sorted_aligned_reads.bam 2> ${fastq_simple_name}_novosort.log

            htseq-count -f bam -o ${fastq_simple_name}_htseq_annote.sam -s ${strandedness} -t gene -i ID ${fastq_simple_name}_sorted_aligned_reads.bam ${params.S288C_R64_annotation_file} 1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log

            sed "s/\t//" ${fastq_simple_name}_htseq_annote.sam > ${fastq_simple_name}_htseq_annote.sam

            samtools view ${fastq_simple_name}_sorted_aligned_reads.bam | paste - b > ${fastq_simple_name}_sorted_aligned_reads_with_annote.bam
            """
        else if (organism == 'H99')
            """
            novoalign -r All -c 8 -o SAM -d ${params.H99_novoalign_index} -f ${fastq_file} 2> ${fastq_simple_name}_novoalign.log | \
            samtools view -bS | \
            novosort - --threads 8 --markDuplicates -o ${fastq_simple_name}_sorted_aligned_reads.bam 2> ${fastq_simple_name}_novosort.log

            htseq-count -f bam -s ${strandedness} -t exon -i gene_id ${fastq_simple_name}_sorted_aligned_reads.bam ${params.H99_annotation_file} 1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log

            sed "s/\t//" ${fastq_simple_name}_htseq_annote.sam > ${fastq_simple_name}_htseq_annote.sam

            samtools view ${fastq_simple_name}_sorted_aligned_reads.bam | paste - b > ${fastq_simple_name}_sorted_aligned_reads_with_annote.bam
            """
  }
