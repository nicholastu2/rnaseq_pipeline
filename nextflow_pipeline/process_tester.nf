#!/usr/bin/env nextflow

process fastqc {
    input:


    script:
    """
    echo "I\'m working!" >> /scratch/mblab/chasem/nextflow_output_tester.txt
    """
}
