#!/usr/bin/env python
import sys
import os
import argparse
from glob import glob
import re

FASTQ_TYPES = ["fastq.gz", "fq.gz"]
FEATURE_TYPE_DICT = {"gff": "gene", "gtf": "CDS"}

def main(argv):

    args = parse_args(argv)
    fastq_path = args.fastq_path
    geno_idx_file = args.genome_index_file
    gene_ann_file = args.gene_annotation_file
    ann_feat_type = args.annotation_feature_type
    strandness = args.strandness
    user_email = args.user_email
    align_only = args.align_only
    run_num = args.run_num
    output_path = os.path.join(args.output_path, 'run_{}'.format(run_num))

    print('...parsing input')
    # Parse default variables
    if output_path is None:
        output_path = fastq_path

    if gene_ann_file is None:
        gene_ann_prefix = ".".join(geno_idx_file.split(".")[:-1])
        gene_ann_file = find_annotation_file(gene_ann_prefix)

    if ann_feat_type is None:
        feat_type = get_feature_type(gene_ann_file)

    if run_num is None:
        run_num = get_run_number(fastq_path)

    print('...writing sbatch job script')
    # Write sbatch script
    fastq_list_file = "job_scripts/fastq_list.txt"
    sbatch_job_file = "job_scripts/mblab_rnaseq.sbatch"
    os.system("mkdir -p sbatch_log/")
    os.system("mkdir -p job_scripts/")

    num_fastqs = write_fastq_list(fastq_path, fastq_list_file)

    write_job_script(sbatch_job_file, output_path,
                     fastq_list_file, num_fastqs,
                     geno_idx_file, gene_ann_file, feat_type,
                     strandness, align_only)

    print('...submitting job')
    # Submit sbatch job
    if user_email is None:
        os.system("sbatch {}".format(sbatch_job_file))
    else:
        os.system("sbatch --mail-type=END,FAIL --mail-user={0} {1}".format(user_email, sbatch_job_file))

    # create a subdirectory in run_#### (where align/counts are stored) for pipeline and annotation version info
    output_subdir_path = os.path.join(output_path, "pipeline_info")
    os.system("mkdir -p {}".format(output_subdir_path))
    # get version info from the module .lua file (see the .lua whatis statements)
    pipeline_info_path = os.path.join(output_subdir_path, 'pipeline_info.txt')
    os.system("module whatis rnaseq_pipeline 2> {}".format(pipeline_info_path))
    with open(pipeline_info_path, "a+") as file:
        file.write("\n")
    os.system("head {} >> pipeline_info_path".format(gene_ann_file))


def parse_args(argv):
    parser = argparse.ArgumentParser(description="This script generates sbatch script and submits sbatch job.")
    parser.add_argument("-f", "--fastq_path", required=True,
                        help="[Required] Directory path of fastq files of type: {}".format(FASTQ_TYPES))
    parser.add_argument("-i", "--genome_index_file", required=True,
                        help="[Required] File path of genome index. This is specific to the aligner.")
    parser.add_argument("-s", "--strandness", required=True,
                        help="[Required] Specify 'yes', 'no', or 'reverse'. For NEB kit, use 'reverse'.")
    parser.add_argument('-o', "--output_path", required=True,
                        help="[Required]  Suggested usage: reports. \n \
                        The topmost directory in which to deposit count, bam and novoalign log files \n \
                        in a automatically generated subdirectory named by the run number.")
    parser.add_argument('-r', '--run_num', required=True,
                        help='[Required] To be used in the event that there is no run number in the fastq_path')
    parser.add_argument("-ga", "--gene_annotation_file", required=True,
                        help="[Required] File path of gene annotation. By default (if not specified), it will look for .gff or .gtf file in the same directory and has same filename as genome index file.")
    parser.add_argument("--annotation_feature_type", default=None,
                        help="[Optional]  Feature type to use for reads counting. By default (if not specified), it will use this dictionary {} based on the annotation file type.".format(
                            FEATURE_TYPE_DICT))
    parser.add_argument("--user_email", default=None,
                        help="[Optional] Email for job status notification.")
    parser.add_argument("--align_only", action="store_true",
                        help="[Optional] Set this flag to only align reads.")
    args = parser.parse_args(argv[1:])
    return args


def find_annotation_file(path_prefix):
    # looks for annotation files
    # Args: the filepath to the files in the annotation dir (/annotation/my_annot.{gtf, gff}
    # Return: the inputted annotation path, if it exists
    # TODO: turn into try/catch error handling
    parsed = None
    for suffix in FEATURE_TYPE_DICT.keys():
        file_path = path_prefix + "." + suffix
        if os.path.exists(file_path):
            parsed = file_path
    if parsed is None:
        sys.exit("Error: annotation file not found.")
    return parsed


def get_feature_type(file_path):
    # get annotation feature type -- this is used if no user input
    # Args: filepath to where .gff and .gtf are
    # Return: feature type of file_path if exists

    suffix = file_path.split(".")[-1]
    if suffix not in FEATURE_TYPE_DICT.keys():
        sys.exit("Error: annotation file has incorrect suffix.")
    return FEATURE_TYPE_DICT[suffix]


def write_fastq_list(dir_path, fastq_list_file):
    # write fastq filepaths in a list stored as a .txt. Used in slurm job script
    # Args: dir_path to fastq_files, out_file to output (where the necessary subdirectories are, in particular job_scripts. see rnaseq pipeline wiki)
    # Returns: the number of fastqs to be aligned/
    # REQUIRES: job_scripts be in directory from which this is called
    write_path = fastq_list_file
    file_paths = []
    for suffix in FASTQ_TYPES:
        file_paths += glob(dir_path + "/*." + suffix)
    with open(write_path, "w") as f:
        for file_path in file_paths:
            f.write("{}\n".format(file_path))
    return len(file_paths)


def get_run_number(fastq_path):
    # extract filepath from fastq filepath. Expecting the directory which stores fastq files to be run_\d*_samples (string on either side does not matter)
    # Args: the cmd line input fastq filepath directory
    # Return: the run number as string (leading zero will be retained)
    try:
        regex = r"run_(\d*)_samples"
        run_num = re.search(regex, fastq_path).group(1)
    except AttributeError:
        print("\nrun_number not found in the fastq_filename input path. \
              Please use the optional cmd line input to enter a run number")
    else:
        return run_num


def write_job_script(job_file, output_path, fastq_list_file, num_fastqs, geno_idx_file, gene_ann_file, feat_type,
                     strandness, align_only):
    # create the job slurm job submission
    # Args: see cmd line input
    # Return: slurm job script

    with open(job_file, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("#SBATCH -N 1\n")
        f.write("#SBATCH --cpus-per-task=8\n")
        f.write("#SBATCH --mem=12G\n")
        f.write("#SBATCH --array=1-{0}%{1}\n".format(num_fastqs, min(num_fastqs, 50)))
        f.write("#SBATCH -D ./\n")
        f.write("#SBATCH -o sbatch_log/mblab_rnaseq_%A_%a.out\n")
        f.write("#SBATCH -e sbatch_log/mblab_rnaseq_%A_%a.err\n")
        f.write("#SBATCH -J mblab_rnaseq\n\n")
        f.write("ml novoalign/3.07.00\n")
        f.write("ml samtools/1.6\n")
        f.write("ml htseq/0.9.1\n")
        f.write("read fastq_file < <( sed -n ${{SLURM_ARRAY_TASK_ID}}p {} ); set -e\n\n".format(fastq_list_file))
        f.write("mkdir -p {}\n".format(output_path))
        f.write(
            "sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; novoalign -c 8 -o SAM -d {0} -f ${{fastq_file}} 2> {1}/${{sample}}_novoalign.log | samtools view -bS > {1}/${{sample}}_aligned_reads.bam\n".format(
                geno_idx_file, output_path))
        f.write(
            "sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; novosort --threads 8 {0}/${{sample}}_aligned_reads.bam > {0}/${{sample}}_sorted_aligned_reads.bam 2> {0}/${{sample}}_novosort.log\n".format(
                output_path))
        if align_only is False:
            if feat_type == 'gene': # this is a bad way of saying "if gff, specify -i ID. Else (it is a gtf) do not specify ID. TODO: This needs to be cleaned up
                f.write(
                    "sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; htseq-count -f bam -i ID -s {1} -t {2} {0}/${{sample}}_sorted_aligned_reads.bam {3} > {0}/${{sample}}_read_count.tsv 2> {0}/${{sample}}_htseq.log\n".format(output_path, strandness, feat_type, gene_ann_file))
            else:
                f.write("sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; htseq-count -f bam -s {1} -t {2} {0}/${{sample}}_sorted_aligned_reads.bam {3} > {0}/${{sample}}_read_count.tsv 2> {0}/${{sample}}_htseq.log\n".format(output_path, strandness, feat_type, gene_ann_file))

if __name__ == "__main__":
    main(sys.argv)