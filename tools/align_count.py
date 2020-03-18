#!/usr/bin/env python
import sys
import os
import argparse
from glob import glob
import re
import datetime
from utils import pathBaseName

FASTQ_TYPES = ["fastq.gz", "fastq", "fq.gz", "fq"]
FEATURE_TYPE_DICT = {"gff": "gene", "gtf": "CDS"}

def main(argv):

    # parse command line input and store as more descriptive variables
    print('...parsing input')
    args = parse_args(argv)
    fastq_path = args.fastq_path
    genome_index_file = args.genome_index_file
    genome_annotation_file = args.gene_annotation_file
    annotation_feature_type = args.annotation_feature_type
    strandness = args.strandness
    user_email = args.user_email
    align_only = args.align_only
    run_num = args.run_num
    output_path = os.path.join(args.output_directory, 'run_{}'.format(run_num))
    if run_num is None:
        run_num = getRunNumber(fastq_path)
    if annotation_feature_type is None:
        annotation_feature_type = getFeatureType(genome_annotation_file)

    # get current datetime
    current_datetime = datetime.datetime.now()
    # get basename of annotation file to name pipeline_info in event that there is more than one organism in the run
    annotation_file_basename = pathBaseName(genome_annotation_file)

    # Create directories to store job scripts and logs, if they dne, and create job script files
    print('...writing sbatch job script')
    # create directories sbatch_log and job_scripts in $PWD if they do not already exist
    os.system("mkdir -p sbatch_log/")
    os.system("mkdir -p job_scripts/")
    # create filenames for job_scripts and sbatch_logs
    fastq_list_file = "job_scripts/run_{}_fastq_list.txt".format(run_num)
    sbatch_job_file = "job_scripts/run_{}_mblab_rnaseq.sbatch".format(run_num)
    # write a list of fastqfiles in the fastq_path to ./job_scripts and store how many as num_fastqs
    num_fastqs = writeFastqList(fastq_path, fastq_list_file)
    # create a slurm submission script and write to ./job_scripts
    writeJobScript(sbatch_job_file, output_path, fastq_list_file, num_fastqs, genome_index_file, genome_annotation_file,
                   annotation_feature_type, strandness, align_only)

    # execute slurm job
    print('...submitting job')
    if user_email is None:
        os.system("sbatch {}".format(sbatch_job_file))
    else:
        os.system("sbatch --mail-type=END,FAIL --mail-user={0} {1}".format(user_email, sbatch_job_file))

    # record information about the pipeline and annotation file in output_path/<annotation_file_basename>_pipeline_info
    output_subdir_path = os.path.join(output_path, "{}_pipeline_info".format(annotation_file_basename))
    os.system("mkdir -p {}".format(output_subdir_path))
    # write version info from the module .lua file (see the .lua whatis statements)
    pipeline_info_path = os.path.join(output_subdir_path, 'pipeline_info.txt')
    os.system("module whatis rnaseq_pipeline 2> {}".format(pipeline_info_path))
    # include the date processed in output_subdir_path/pipeline_into.txt
    with open(pipeline_info_path, "a+") as file:
        file.write("\n")
        file.write('Date processed: {:%Y-%m-%d %H:%M:%S}'.format(current_datetime))
        file.write("\n")
    # include the head of the gff/gtf, also
    os.system("head {} >> {}".format(annotation_feature_type, pipeline_info_path))

def parse_args(argv):
    parser = argparse.ArgumentParser(description="This script generates sbatch script and submits sbatch job.")
    parser.add_argument("-f", "--fastq_path", required=True,
                        help="[Required] Directory path of fastq files of type: {}".format(FASTQ_TYPES))
    parser.add_argument("-i", "--genome_index_file", required=True,
                        help="[Required] File path of genome index. This is specific to the aligner.")
    parser.add_argument("-ga", "--gene_annotation_file", required=True,
                        help="[Required] File path of gene annotation. By default (if not specified), it will look for \
                        .gff or .gtf file in the same directory and has same filename as genome index file.")
    parser.add_argument("-s", "--strandness", required=True,
                        help="[Required] Specify 'yes', 'no', or 'reverse'. For NEB kit, use 'reverse'.")
    parser.add_argument('-o', "--output_directory", required=True,
                        help="[Required]  Suggested usage: reports. \n \
                        The topmost directory in which to deposit count, bam and novoalign log files \n \
                        in a automatically generated subdirectory named by the run number.")
    parser.add_argument('-rn', '--run_num',
                        help='[Optional] To be used in the event that there is no run number in the fastq_path')
    parser.add_argument("--user_email", default=None,
                        help="[Optional] Email for job status notification.")
    parser.add_argument("--annotation_feature_type", default=None,
                        help="[Optional]  Feature type to use for reads counting. By default (if not specified), it will \
                        use this dictionary {} based on the annotation file type.".format(FEATURE_TYPE_DICT))
    parser.add_argument("--align_only", action="store_true",
                        help="[Optional] Set this flag to only align reads.")

    args = parser.parse_args(argv[1:])
    return args

def getRunNumber(fastq_path):
    """
    extract run number from -f input. this *should be* to a file called sequence/run_####_samples
    :param fastq_path
    :returns: the run number as string (leading zero will be retained)
    """
    try:
        regex = r"run_(\d*)_samples"
        run_num = re.search(regex, fastq_path).group(1)
    except AttributeError:
        print("\nrun_number not found in the fastq_filename input path. \
              Please use the optional cmd line input to enter a run number")
    else:
        return run_num

def getFeatureType(file_path):
    """
    :param file_path: filepath to .gff or .gtf
    :returns: feature type according to FEATURE_TYPE_DICT
    """
    suffix = file_path.split(".")[-1]
    if suffix not in FEATURE_TYPE_DICT.keys():
        sys.exit("Error: annotation file has incorrect suffix.")
    return FEATURE_TYPE_DICT[suffix]

def writeFastqList(dir_path, fastq_list_file):
    """
    write fastq filepaths in a list stored as a .txt. Used in slurm job script
    :param dir_path:
    :param fastq_list_file:
    :returns: the number of fastqs to be aligned/
    """

    write_path = fastq_list_file
    file_paths = []
    for suffix in FASTQ_TYPES:
        file_paths += glob(dir_path + "/*." + suffix)
    with open(write_path, "w") as f:
        for file_path in file_paths:
            f.write("{}\n".format(file_path))
    return len(file_paths)

def writeJobScript(job_file, output_path, fastq_list_file, num_fastqs, genome_index_file, genome_annotation_file, feature_type,
                   strandness, align_only):
    """
    Write slurm job script to job_file (which is $PWD/job_scripts
    :param job_file: path to $PWD/job_scripts (see main method)
    :param output_path: path to output_dir (see main method)
    :param fastq_list_file: path to job_scripts/<run_number>_fastq_list.txt
    :param num_fastqs: number of lines in fastq_list_file
    :param genome_index_file: path to index file created from genome of given organism by novoalign
    :param genome_annotation_file: path to genome annotation file, either .gff or .gtf
    :param feature_type: feature type extracted based on FEATURE_TYPE_DICT (see top of script)
    :param strandness: cmd line input from user regarding whether library prep is stranded
    :param align_only: boolean flag allowing cmd line input for alignment only, no htseq
    :returns: None
    """

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
                genome_index_file, output_path))
        f.write(
            "sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; novosort --threads 8 {0}/${{sample}}_aligned_reads.bam > {0}/${{sample}}_sorted_aligned_reads.bam 2> {0}/${{sample}}_novosort.log\n".format(
                output_path))
        if align_only is False:
            if feature_type == 'gene': # this is a bad way of saying "if gff, specify -i ID. Else (it is a gtf) do not specify ID. TODO: This needs to be cleaned up
                f.write(
                    "sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; htseq-count -f bam -i ID -s {1} -t {2} {0}/${{sample}}_sorted_aligned_reads.bam {3} > {0}/${{sample}}_read_count.tsv 2> {0}/${{sample}}_htseq.log\n".format(output_path, strandness, feature_type, genome_annotation_file))
            else:
                f.write("sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; htseq-count -f bam -s {1} -t {2} {0}/${{sample}}_sorted_aligned_reads.bam {3} > {0}/${{sample}}_read_count.tsv 2> {0}/${{sample}}_htseq.log\n".format(output_path, strandness, feature_type, genome_annotation_file))


if __name__ == "__main__":
    main(sys.argv)
