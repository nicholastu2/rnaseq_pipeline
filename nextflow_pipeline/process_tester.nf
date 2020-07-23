// need to figure out pipeline info process

process htseq_count {

    scratch true
    executor "slurm"
    cpus 8
    memory "20G"
    beforeScript "ml samtools htseq"
    stageInMode "copy"
    stageOutMode "move"
    publishDir "$params.align_count_results/$run_directory/logs", mode:"copy", overwite: true, pattern: "*.log"
    publishDir "$params.align_count_results/$run_directory/count", mode:"copy", overwite: true, pattern: "*_read_count.tsv"
    publishDir "$params.align_count_results/$run_directory/align", mode:"copy", overwite: true, pattern: "*_sorted_aligned_reads_with_annote.bam"


    input:
        tuple val(run_directory), val(fastq_simple_name), val(organism), val(strandedness), file(sorted_bam) from bam_align_ch
    output:
        tuple val(run_directory), val(fastq_simple_name), val(organism), val(strandedness), file("${fastq_simple_name}_sorted_aligned_reads_with_annote.bam") into bam_align_with_htseq_annote_ch
        tuple val(run_directory), val(fastq_simple_name), file("${fastq_simple_name}_read_count.tsv") into htseq_count_ch
        tuple val(run_directory), val(fastq_simple_name), file("${fastq_simple_name}_htseq.log") into htseq_log_ch
        tuple val(run_directory), val(organism), val(strandedness) pipeline_info_ch

    script:
        if (organism == 'S288C_R64')
            """
            htseq-count -f bam \\
                        -o ${fastq_simple_name}_htseq_annote.sam \\
                        -s ${strandedness} \\
                        -t gene \\
                        -i ID \\
                        ${sorted_bam} \\
                        ${params.S288C_R64_annotation_file} \\
                        1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log

            sed "s/\t//" ${fastq_simple_name}_htseq_annote.sam > ${fastq_simple_name}_no_tab_sam.sam

            samtools view --threads 8 ${sorted_bam} | \\
            paste - ${fastq_simple_name}_no_tab_sam.sam | \\
            samtools view --threads 8 -bS -T ${params.S288C_R64_genome} > ${fastq_simple_name}_sorted_aligned_reads_with_annote.bam

            """
        else if (organism == 'KN99' && strandedness == 'reverse')
            """
            htseq-count -f bam \\
                        -o ${fastq_simple_name}_htseq_annote.sam \\
                        -s ${strandedness} \\
                        -t exon \\
                        -i gene \\
                        ${sorted_bam} \\
                        ${params.KN99_annotation_file} \\
                        1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log

            sed "s/\t//" ${fastq_simple_name}_htseq_annote.sam > ${fastq_simple_name}_no_tab_sam.sam

            samtools view --threads 8 ${sorted_bam} | \\
            paste - ${fastq_simple_name}_no_tab_sam.sam | \\
            samtools view --threads 8 -bS -T ${params.KN99_genome} > ${fastq_simple_name}_sorted_aligned_reads_with_annote.bam

            """
        else if (organism == 'KN99' && strandedness == 'no')
            """
            htseq-count -f bam \\
                        -o ${fastq_simple_name}_htseq_annote.sam \\
                        -s ${strandedness} \\
                        -t exon \\
                        -i gene \\
                        ${sorted_bam} \\
                        ${params.KN99_annotation_file_no_strand} \\
                        1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log

            sed "s/\t//" ${fastq_simple_name}_htseq_annote.sam > ${fastq_simple_name}_no_tab_sam.sam

            samtools view --threads 8 ${sorted_bam} | \\
            paste - ${fastq_simple_name}_no_tab_sam.sam | \\
            samtools view --threads 8 -bS -T ${params.KN99_genome} > ${fastq_simple_name}_sorted_aligned_reads_with_annote.bam

            """
        else if (organism == 'H99')
            """
            htseq-count -f bam \\
                        -s ${strandedness} \\
                        -t exon \\
                        ${sorted_bam} \\
                        ${params.H99_annotation_file} \\
                        1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log

            sed "s/\t//" ${fastq_simple_name}_htseq_annote.sam > ${fastq_simple_name}_no_tab_sam.sam

            samtools view --threads 8 ${sorted_bam} | \\
            paste - ${fastq_simple_name}_no_tab_sam.sam | \\
            samtools view --threads 8 -bS -T ${params.H99_genome} > ${fastq_simple_name}_sorted_aligned_reads_with_annote.bam

            """
}

process writePipelineInfo {

    executor "local"
    beforeScript "ml rnaseq_pipeline"

    input:
        tuple val(run_directory), val(organism), val(strandedness) from pipeline_info_ch

    script:
        """
        #!/usr/bin/env python

        from rnaseq_tools.OrganismDataObject import OrganismDataObject
        from rnaseq_tools import utils

        # instantiate OrganismDataObject (see brentlab rnaseq_pipeline)
        od = OrganismDataObject(organism = ${organism}, interactive=True)

        # create pipeline_info subdir of in rnaseq_pipeline/align_count_results/${organism}_pipeline_info
        pipeline_info_subdir_path = os.path.join(od.align_count_results, ${run_directory}, "${organism}_pipeline_info")
        utils.mkdirp(pipeline_info_subdir_path)

        # write version info from the module .lua file (see the .lua whatis statements)
        pipeline_info_txt_file_path = os.path.join(pipeline_info_subdir_path, 'pipeline_info.txt')
        cmd_pipeline_info = "module whatis rnaseq_pipeline 2> {}".format(pipeline_info_txt_file_path)
        utils.executeSubProcess(cmd_pipeline_info)

        # include the date processed in pipeline_info_subdir_path/pipeline_into.txt
        with open(pipeline_info_txt_file_path, "a+") as file:
            file.write("\n")
            current_datetime = od.year_month_day + '_' + utils.hourMinuteSecond()
            file.write('Date processed: %s' % current_datetime)
            file.write("\n")

        # set annotation_file
        if ${organism} == 'KN99':
            if ${stranded} == 'no':
                annotation_file = od.annotation_file_no_strand
        else:
            annotation_file = od.annotation_file
        # include the head of the gff/gtf in pipeline_info
        cmd_annotation_info = "head {} >> {}".format(annotation_file, pipeline_info_txt_file_path)
        utils.executeSubProcess(cmd_annotation_info)

        # TODO: try copying nextflow jobscript to pipeline_info

        """
}
