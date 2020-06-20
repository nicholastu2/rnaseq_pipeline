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
        tuple val(run_directory), val(fastq_simple_name), val(organism), val(strandedness), file("${fastq_simple_name}_aligned_reads.sam"), file("${fastq_simple_name}_novoalign.log") into sam_align_ch

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

processs splitAlignmentFilesByMappingType {
  scratch true
  executor "slurm"
  memory "12G"
  beforeScript "ml samtools"
  stageInMode "copy"
  stageOutMode "move"
  afterScript "rm ${alignment_sam}" // figure out how to actually delete these


}

// blast unmapped reads
process blastUnmapped {
// convert bam to fasta, blast, publish results immediately
// use < bamToFasta into blast > blast_results.out

}

// combine this step with novosort
// update to handle split output
process convertSamToBam {
  executor "slurm"
  memory "12G"
  beforeScript "ml samtools"
  stageInMode "copy"
  stageOutMode "move"
  afterScript "rm ${alignment_sam}" // fix this to handle split files

  input:
    tuple val(run_directory), val(fastq_simple_name), val(organism), val(strandedness), file(alignment_sam), file(novoalign_log) from sam_align_ch

  output:
    tuple val(run_directory), val(fastq_simple_name), val(organism), val(strandedness), file("${fastq_simple_name}_aligned_reads.bam"), file(novoalign_log) into bam_align_ch

  script:
     """
     samtools view -bS ${alignment_sam}> ${fastq_simple_name}_aligned_reads.bam
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
  //publishDir --> just publish novosort_log now

    input:
      tuple val(run_directory), val(fastq_simple_name), val(organism), val(strandedness), file(alignment_sam), file(novoalign_log) from sam_align_ch
    output:
      tuple val(run_directory), val(fastq_simple_name), val(organism), val(strandedness), file(novoalign_log), file("${fastq_simple_name}_sorted_aligned_reads.bam"), file("${fastq_simple_name}_novosort.log") into sorted_bam_align_ch

    // check what happens with the --index option
    // see http://www.novocraft.com/documentation/novosort-2/ for novosort options

    script:
      """
      samtools view -bS ${alignment_sam} | novosort --threads 8 --markDuplicates --index --output ${fastq_simple_name}_sorted_aligned_reads.bam 2> ${fastq_simple_name}_novosort.log
      """
}

// needs to be updated to output the --out bam info
process htseqCount {
  scratch true
  executor "slurm"
  memory "12G"
  stageInMode "copy"
  stageOutMode "move"
  beforeScript "ml htseq"
  publishDir "$params.align_count_results/$run_directory", mode:"move", overwite: true, pattern: "*_novoalign.log"
  publishDir "$params.align_count_results/$run_directory", mode:"move", overwite: true, pattern: "*_novosort.log"
  publishDir "$params.align_count_results/$run_directory", mode:"move", overwite: true, pattern: "*_sorted_aligned_reads.bam"
  publishDir "$params.align_count_results/$run_directory", mode:"move", overwite: true, pattern: "*_read_count.tsv"
  publishDir "$params.align_count_results/$run_directory", mode:"move", overwite: true, pattern: "*_htseq.log"

    input:
      tuple val(run_directory), val(fastq_simple_name), val(organism), val(strandedness), file(novoalign_log), file(sorted_alignment_bam), file(novosort_log) from sorted_bam_align_ch
    output:
      tuple file(novoalign_log), file(novosort_log), file(sorted_alignment_bam), file("${fastq_simple_name}_read_count.tsv"), file("${fastq_simple_name}_htseq.log") into align_count_output_ch
      tuple val(run_directory), val(organism) into pipeline_info_ch

    script:
    // note: this is set for gff3 format. eg: (tab delimited)
    // CP022321.1	EuPathDB	exon	1767728	1768158	.	+	.	ID=exon_CKF44_00681-E1;Parent=CKF44_00681-t42_2,CKF44_00681-t42_1
    // hence -t exon -i ID
      if (organism == 'KN99')
        """
        htseq-count -f bam -s ${strandedness} -t exon -i Parent --additional-attr ID ${sorted_alignment_bam} ${params.KN99_annotation_file} 1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log
        """
      else if (organism == 'S288C_R64')
        """
        htseq-count -f bam -s ${strandedness} -t gene -i ID ${sorted_alignment_bam} ${params.S288C_R64_annotation_file} 1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log
        """
      else if (organism == 'H99')
        """
        htseq-count -f bam -s ${strandedness} -t exon -i gene_id ${sorted_alignment_bam} ${params.H99_annotation_file} 1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log
        """
      else
        error "Invalid Organism"

}

process featureToBam{
//input tracer data + sorted bam
// output tracer data + sorted_bam_plus_features

}

// collapse parent to gene eg CKF44_00681
process KN99TranscriptToGene {

}

process qualAssess1 {

  scratch true
  executor "slurm"
  memory "12G"
  stageInMode "copy"
  stageOutMode "move"
  executor "slurm"
  beforeScript "ml rnaseq_pipeline"
  // write everything out

  input:
     // collect logs, alignment, htseq output gather by run runNumber from input_ch.collect().groupTuple()

  output:
      // everything + quality assess

  script:
      """
      #!/usr/bin/env python

      import sys
      import os
      import argparse
      from rnaseq_tools.QualityAssessmentObject import QualityAssessmentObject
      from rnaseq_tools import utils

      # for ordering columns below. genotype_1_coverage and genotype_2_coverage added if coverage_check is passed
      column_order = ['LIBRARY_SIZE', 'TOTAL_ALIGNMENT', 'UNIQUE_ALIGNMENT', 'MULTI_MAP', 'NO_MAP', 'HOMOPOLY_FILTER',
                      'READ_LENGTH_FILTER', 'NOT_ALIGNED_TOTAL', 'WITH_FEATURE', 'NO_FEATURE', 'FEATURE_ALIGN_NOT_UNIQUE',
                      'AMBIGUOUS_FEATURE', 'TOO_LOW_AQUAL', 'GENOTYPE_1_COVERAGE', 'GENOTYPE_2_COVERAGE']

      # create QualityAssessmentObject (see rnaseq_tools in rnaseq_pipeline repo)
      qa = QualityAssessmentObject(nextflow_list_of_files=${bam_files} + ${log_files} + {htseq_counts},
                                   log_suffix_list=["_novoalign.log", "_read_count.tsv"],
                                   coverage_check_flag=True,
                                   query_path = ${query_path},
                                   interactive=True)
                                   print('...compiling alignment information')

      # create dataframes storing the relevant alignment and count metadata from the novoalign and htseq logs
      qual_assess_1_df = qa.compileData()

      # re_order columns
      qual_assess_1_df = qual_assess_1_df[column_order]

      # write out (name can be generic -- stored in unique work directory)
      qual_assess_1_df.to_csv('quality_assess_1.csv', index_label = 'FASTQFILENAME')

      """
}


process failedPerturbedGeneBrowserShot{



}
