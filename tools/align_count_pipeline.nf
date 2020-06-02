#!/usr/bin/env nextflow

// see http://nextflow-io.github.io/patterns/index.html

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
        tuple val(containing_directory), file(fastq_list) from samples_channel

    script:
    """
    echo your_command --sample $sampleId --reads $read1 $read2
    """
}

process multi_qc {
    input:
    set sampleId, file(read1), file(read2) from samples_channel

    script:
    """
    echo your_command --sample $sampleId --reads $read1 $read2
    """
}

process qualimap {
    input:
    set sampleId, file(read1), file(read2) from samples_channel

    script:
    """
    echo your_command --sample $sampleId --reads $read1 $read2
    """
}


process lts_to_scratch {
    cache = 'false'
    executor = 'local'
    beforeScript = 'ml rnaseq_pipeline'

    input:
        set sampleId, file(read1), file(read2) from samples_channel

    output:
        tuple 'scratch_path',
        file param.scratch_sequence/${}.fastq.gz into fastq_process_list

    script:
    """
    // extract parent directory from source, make directory of same name to deposit file in destination
    rsync_copy_nextflow.py -s --sample $fastq_file -d  param.scratch_sequence //somehow this needs to go to a run_#### file in scratch sequence b/c fastqfilenames may not be unique
    """
}

process novoalign {

    beforeScript = 'ml novoalign'
    afterScript // remove fastq here!
    cpus = 8
    executor = 'slurm'

    input:
    set sampleId, file(read1), file(read2) from samples_ch

    script:
    if ($organism == 'KN99') //figure out how to redirect output
    """
    novoalign -r All -c 8 -o SAM -d param.KN99_genome_index -f ${{fastq_file}} 2> {1}/${{sample}}_novoalign.log".format(output_path))
    """
    else if ($organism == 'S288C_R64')
    """
    novoalign -r All -c 8 -o SAM -d param.S288C_R64_genome_index -f ${{fastq_file}} 2> {1}/${{sample}}_novoalign.log".format(output_path))
    """
    else if ($oganism == 'H99')
    """
    novoalign -r All -c 8 -o SAM -d param.H99_genome_index -f ${{fastq_file}} 2> {1}/${{sample}}_novoalign.log".format(output_path))
    """
    else
    error "Invalid Organism"
}

process convert sam_to_bam {
    // output from novoalign shunted into novosort

    script:
    """
    samtools view -bS > {1}/${{sample}}_aligned_reads.bam\n".format(output_path)
    """

 }

process novosort {

    publishDir
    beforeScript = 'ml novoalign'
    afterScript // remove the align file
    cpus = 8
    executor = 'slurm'

    input:
    set sampleId, file(read1), file(read2) from samples_ch

    script:
    """
    novosort --threads 8 {0}/${{sample}}_aligned_reads.bam > {0}/${{sample}}_sorted_aligned_reads.bam 2> {0}/${{sample}}_novosort.log\n".format(output_path))
    """
}

process htseq_count {

    publishDir
    beforeScript = 'ml htseq'
    input:
    set sampleId, file(read1), file(read2) from samples_ch

    script:
    if ($organism == 'KN99')
    """
    htseq-count -f bam -s $strandedness -t exon ${{sample}} param.KN99_annotation_file > {0}/${{sample}}_read_count.tsv 2> {0}/${{sample}}_htseq.log\n".format(output_path)
    """
    else if ($organism == 'S288C_R64')
    """
    htseq-count -f bam -s $strandedness -t gene ${{sample}} param.S288C_R64_annotation_file > {0}/${{sample}}_read_count.tsv 2> {0}/${{sample}}_htseq.log\n".format(output_path)
    """
    else if ($oganism == 'H99')
    """
    htseq-count -f bam -s $strandedness -t exon ${{sample}} param.H99_annotation_file > {0}/${{sample}}_read_count.tsv 2> {0}/${{sample}}_htseq.log\n".format(output_path)
    """
    else
    error "Invalid Organism"
}

process perturbed_gene_coverage {

    publishDir

    input:
    set sampleId, file(read1), file(read2) from samples_ch

    script:
    """
    echo your_command --sample $sampleId --reads $read1 $read2
    """
}

process scratch_to_lts {
    input: --> multiple inputs
    set sampleId, file(read1), file(read2) from samples_ch
    afterScript // remove the files from scratch space


    script:
    """
    echo your_command --sample $sampleId --reads $read1 $read2
    """
}


params.in = "$baseDir/data/sample.fa"

/*
 * split a fasta file in multiple files
 */
process splitSequences {

    input:
    path 'input.fa' from params.in

    output:
    path 'seq_*' into records

    conditional

    """
    awk '/^>/{f="seq_"++d} {print > f}' < input.fa
    """
}

/*
 * Simple reverse the sequences
 */
process reverse {

    input:
    path x from records

    output:
    stdout into result

    """
    cat $x | rev
    """
}

/*
 * print the channel content
 */
result.subscribe { println it }