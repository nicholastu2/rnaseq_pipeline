#!/usr/bin/env python
"""
    A script to run a standard align or align+count job. If you wish to use a non standard genome, feature, etc.
    Then you may find it easier to do so using rnaseq_tools/OrganismDataObject + rnaseq_tools/SbatchWriterObject.
    The functions in rnaseq_tools/utils may also be of use (especially getFileListFromDirectory)

    for help:
        align_count.py -h

    usage: align_count.py -f /path/to/directory/in/scratch/run_####_samples -g <an organism in genome_files> -s yes/no/reverse -o output_dir --user_email youremail@mail.com

    written by: yiming kang https://github.com/yiming-kang/rnaseq_pipeline
    modified by: chase mateusiak chasem@wustl.edu chase.mateusiak@gmail.com
    last update: 05/09/2020
"""

import sys
import os
import argparse
from rnaseq_tools.OrganismDataObject import OrganismData
from rnaseq_tools.SbatchWriterObject import SbatchWriter
from rnaseq_tools import utils

# TODO: re-break up main script into functions for readability

def main(argv):
    # parse command line input and store as more descriptive variables
    print('...parsing input')
    args = parse_args(argv)
    try:
        if not os.path.isdir(args.fastq_path):
            raise NotADirectoryError('FastqDirectoryDoesNotExist')
    except NotADirectoryError:
        print('The path to %s for the raw fastq_files does not exist. Correct and re-submit. Remember this directory cannot be in long term storage')
    print('...creating OrganismDataObject')
    od = OrganismData(organism=args.organism,
                      fastq_path=args.fastq_path,
                      strandness=args.strandness,
                      email=args.user_email,
                      run_number=utils.getRunNumber(args.fastq_path))
    # check directory structure and set organism data (see OrganismData.setOrganismData())
    od.setOrganismData()
    # create logger for this script if od logger is set
    if os.path.isfile(od.log_file_path):
        logger = utils.createLogger(od.log_file_path, 'align_count.py')
    else:
        logger = utils.createStdOutLogger(name='align_count_logger')

    # add attribute output_dir
    od.output_dir = os.path.join(args.output_directory, 'run_{}'.format(od.run_number))
    # store align_only flag from cmd line
    align_only = args.align_only

    print('...extracting list of fastq files to process')
    fastq_list_file = "{}/run_{}_fastq_list.txt".format(od.job_scripts, od.run_number)
    logger.info('The fastq list file path is %s' %fastq_list_file)
    # extract all files with the extensions in the list from od.fastq_path
    fastq_file_list = utils.getFileListFromDirectory(od.fastq_path, ["fastq.gz", "fastq", "fq.gz", "fq"])
    # store length of list
    num_fastqs = len(fastq_file_list)
    # write list to file
    with open(fastq_list_file) as file:
        for fastq_file_basename in fastq_file_list:
            file.write("{}\n".format(fastq_file_basename))
    if not os.path.isfile(fastq_list_file):
        sys.exit("list of fastq files at %s does not exist" %fastq_list_file)
    else:
        print('list of fastq files may be found at %s' %fastq_list_file)

    print('...writing sbatch job_script')
    # create path for sbatch job_script
    sbatch_job_script_path = "{}/run_{}_mblab_rnaseq.sbatch".format(od.job_scripts, od.run_number)
    logger.info('sbatch job script path is %s' %sbatch_job_script_path)
    # create a slurm submission script and write to ./job_scripts
    SbatchWriter.writeAlignCountJobScript(sbatch_job_script_path, od.output_dir, fastq_list_file, num_fastqs, od.novoalign_index,
                                          od.annotation_file, od.feature_type, od.strandness, align_only)
    if not os.path.isfile(sbatch_job_script_path):
        sys.exit('sbatch job_script does not exist at path %s' %sbatch_job_script_path)
    else:
        print('sbatch script may be found at %s' %sbatch_job_script_path)

    # submit sbatch job
    print('...submitting sbatch job')
    if od.email is None:
        cmd = "sbatch %s" %sbatch_job_script_path
        utils.executeSubProcess(cmd)
    else:
        cmd = "sbatch --mail-type=END,FAIL --mail-user=%s %s" %(od.email, sbatch_job_script_path)
        utils.executeSubProcess(cmd)

    print('\nannotation and pipeline information recorded in {}/run_{}/{}'.format(od.output_dir, od.run_number,
                                                                                  'pipeline_info'))
    output_subdir_path = os.path.join(od.output_dir, "{}_pipeline_info".format(od.organism))
    utils.mkdirp(output_subdir_path)

    # write version info from the module .lua file (see the .lua whatis statements)
    pipeline_info_path = os.path.join(output_subdir_path, 'pipeline_info.txt')
    cmd_pipeline_info = "module whatis rnaseq_pipeline 2> {}".format(pipeline_info_path)
    utils.executeSubProcess(cmd_pipeline_info)
    # include the date processed in output_subdir_path/pipeline_into.txt
    with open(pipeline_info_path, "a+") as file:
        file.write("\n")
        current_datetime = od.year_month_day + '_' + utils.hourMinuteSecond()
        file.write('Date processed: {:%Y-%m-%d %H:%M:%S}'.format(current_datetime))
        file.write("\n")

    # include the head of the gff/gtf, also
    cmd_annotation_info = "head {} >> {}".format(od.annotation_file, pipeline_info_path)
    utils.executeSubProcess(cmd_annotation_info)


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Generate and submit an sbatch script to align-only or align+count a directory of fastq_files")
    parser.add_argument("-f", "--fastq_path", required=True,
                        help="[REQUIRED] Directory path of fastq files.\n"
                             "Must be stored in a directory that begins with run_####_'n"
                             "All fastq files must have extension .fastq.gz")
    parser.add_argument("-g", "--organism", required=True,
                        help="[REQUIRED] The organism. This corresponds to the organisms in genome_files. Currently\n"
                             "The options are H99, KN99 and S288C_R64.")
    parser.add_argument("-s", "--strandness", required=True,
                        help="[REQUIRED] Specify 'yes', 'no', or 'reverse'. For NEB kit, use 'reverse'.")
    parser.add_argument('-o', "--output_directory", required=True,
                        help="[REQUIRED]  Suggested usage: reports.\n"
                             "The topmost directory in which to deposit a directory called \'run_####\'\n"
                             "which will store count, alignment and log files")
    parser.add_argument("--user_email", default=None,
                        help="[OPTIONAL] Email for job status notification.")
    parser.add_argument("--align_only", action="store_true",
                        help="[OPTIONAL] Set this flag to only align reads.")

    args = parser.parse_args(argv[1:])
    return args


if __name__ == "__main__":
    main(sys.argv)
