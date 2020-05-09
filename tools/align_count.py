#!/usr/bin/env python
import sys
import os
import argparse
from glob import glob
import datetime
from rnaseq_tools.OrganismDataObject import OrganismData
from rnaseq_tools import utils

FASTQ_TYPES = ["fastq.gz", "fastq", "fq.gz", "fq"]


def main(argv):
    # parse command line input and store as more descriptive variables
    print('...parsing input')
    args = parse_args(argv)
    # create OrganismData object. This will simultaneously check to ensure that the user has the required file structure
    # in their scratch space in addition to containing the organism attributes (genome, etc)
    sd = OrganismData(organism=args.organism,
                      fastq_path=args.fastq_path,
                      strandness=args.strandness,
                      email=args.user_email,
                      run_number=utils.getRunNumber(args.fastq_path))

    # ensure that standard directory structure is set
    sd.standardDirectoryStructure()
    # create StandardData logger
    sd.setStandardDataLogger()
    sd.output_dir = os.path.join(args.output_directory, 'run_{}'.format(sd.run_number))
    # store align_only flag from cmd line
    align_only = args.align_only
    # get current datetime
    current_datetime = datetime.datetime.now()

    print('...writing sbatch job script')
    # create filenames for job_scripts and sbatch_logs
    fastq_list_file = "{}/run_{}_fastq_list.txt".format(sd.job_scripts, sd.run_number)
    sbatch_job_file = "{}/run_{}_mblab_rnaseq.sbatch".format(sd.job_scripts, sd.run_number)
    # write a list of fastqfiles in the fastq_path to ./job_scripts and store how many as num_fastqs
    num_fastqs = writeFastqList(sd.fastq_path, fastq_list_file)
    # create a slurm submission script and write to ./job_scripts
    writeJobScript(sbatch_job_file, sd.output_dir, fastq_list_file, num_fastqs, sd.novoalign_index,
                   sd.annotation_file, sd.feature_type, sd.strandness, align_only)

    # execute slurm job
    print('...submitting job')
    if sd.email is None:
        os.system("sbatch {}".format(sbatch_job_file))
    else:
        os.system("sbatch --mail-type=END,FAIL --mail-user={0} {1}".format(sd.email, sbatch_job_file))

    print('\nannotation and pipeline information recorded in {}/run_{}/{}'.format(sd.output_dir, sd.run_number,
                                                                                  'pipeline_info'))
    output_subdir_path = os.path.join(sd.output_dir, "{}_pipeline_info".format(sd.organism))
    utils.mkdirp(output_subdir_path)
    # write version info from the module .lua file (see the .lua whatis statements)
    pipeline_info_path = os.path.join(output_subdir_path, 'pipeline_info.txt')
    cmd_pipeline_info = "module whatis rnaseq_pipeline 2> {}".format(pipeline_info_path)
    utils.executeSubProcess(cmd_pipeline_info)
    # include the date processed in output_subdir_path/pipeline_into.txt
    with open(pipeline_info_path, "a+") as file:
        file.write("\n")
        file.write('Date processed: {:%Y-%m-%d %H:%M:%S}'.format(current_datetime))
        file.write("\n")
    # include the head of the gff/gtf, also
    cmd_annotation_info = "head {} >> {}".format(sd.annotation_file, pipeline_info_path)
    utils.executeSubProcess(cmd_annotation_info)


def parse_args(argv):
    parser = argparse.ArgumentParser(description="This script generates sbatch script and submits sbatch job.")
    parser.add_argument("-f", "--fastq_path", required=True,
                        help="[Required] Directory path of fastq files.\n"
                             "Must be stored in a directory that begins with run_####_'n"
                             "All fastq files must have extension .fastq.gz")
    parser.add_argument("-g", "--organism", required=True,
                        help="[REQUIRED] The organism. This corresponds to the organisms in genome_files. Currently\n"
                             "The options are H99, KN99 and S288C_R64.")
    parser.add_argument("-s", "--strandness", required=True,
                        help="[Required] Specify 'yes', 'no', or 'reverse'. For NEB kit, use 'reverse'.")
    parser.add_argument('-o', "--output_directory", required=True,
                        help="[Required]  Suggested usage: reports.\n"
                             "The topmost directory in which to deposit a directory called \'run_####\'\n"
                             "which will store count, alignment and log files")
    parser.add_argument("--user_email", default=None,
                        help="[Optional] Email for job status notification.")
    parser.add_argument("--align_only", action="store_true",
                        help="[Optional] Set this flag to only align reads.")

    args = parser.parse_args(argv[1:])
    return args


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

# TODO: test and then incorporate SbatchWriterObject
def writeJobScript(job_file, output_path, fastq_list_file, num_fastqs, genome_index_file, genome_annotation_file,
                   feature_type,
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

        f.write("sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; novoalign -c 8 -o SAM -d {0} "
                "-f ${{fastq_file}} 2> {1}/${{sample}}_novoalign.log | samtools view -bS > "
                "{1}/${{sample}}_aligned_reads.bam\n".format(genome_index_file, output_path))

        f.write("sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; novosort --threads 8 "
                "{0}/${{sample}}_aligned_reads.bam > {0}/${{sample}}_sorted_aligned_reads.bam 2> "
                "{0}/${{sample}}_novosort.log\n".format(output_path))

        if not align_only:
            if feature_type == 'gene':  # this is a messy way of saying "if gff, specify -i ID. TODO: This needs to be cleaned up
                f.write("sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; "
                        "htseq-count -f bam -i ID -s {1} -t {2} {0}/${{sample}}_sorted_aligned_reads.bam "
                        "{3} > {0}/${{sample}}_read_count.tsv 2> {0}/${{sample}}_htseq.log\n".format(output_path,
                                                                                                     strandness,
                                                                                                     feature_type,
                                                                                                     genome_annotation_file))
            else:  # else (it is a gtf) do not specify ID
                f.write("sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; htseq-count -f bam -s {1} "
                        "-t {2} {0}/${{sample}}_sorted_aligned_reads.bam {3} > "
                        "{0}/${{sample}}_read_count.tsv 2> {0}/${{sample}}_htseq.log\n".format(output_path, strandness,
                                                                                               feature_type,
                                                                                               genome_annotation_file))


if __name__ == "__main__":
    main(sys.argv)
