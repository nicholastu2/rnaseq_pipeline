#!/usr/bin/env python

"""
Purpose: copy count files from /lts into /scratch/mblab/$USER for analysis
Usage: create_experiment.py -qs <your_query>.csv -o . -n <exp_name>
"""
import os
import sys
import pandas as pd
import argparse
from shutil import copy2 as cp
from utils import *

COUNT_LTS = '/lts/mblab/Crypto/rnaseq_data/align_expr'

def main(argv):

    # store suffixes of the files we wish to move
    count_suffix = '_read_count.tsv'
    novoalign_log_suffix = '_novoalign.log'

    # parse cmd line arguments
    args = parseArgs(argv)

    # read in query sheet (The path to the result of a query against the metadata base using queryDB)
    query = pd.read_csv(args.query_sheet)

    # create a directory for the experiment
    dest_dir = createExperimentDir(args.output_location, args.experiment_name)

    # get list of count files
    count_file_list = filepathList(query, count_suffix)
    # get list of novoalign logs
    novoalign_log_list = filepathList(query, novoalign_log_suffix)
    # concat the lists together
    file_list = count_file_list + novoalign_log_list

    # move the files from /lts to the output directory (generally the user's scratch)
    moveFiles(file_list, dest_dir, len(query))


def parseArgs(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-qs', '--query_sheet', required=True,
                        help='query the database using queryDB.py for the libraries you wish to analyze')
    parser.add_argument('-o', '--output_location', required=True,
                        help='Suggested usage: /scratch/$USER or /scratch/$USER/rnaseq  The location where you want \
                            to put the directory containing the counts for your analysis')
    parser.add_argument('-n', '--experiment_name', required=True,
                        help='The name of this experiment. This will be used as the subdirectory of output')
    return parser.parse_args(argv[1:])

def createExperimentDir(output, exp_name):
    # create experiment directory to store copied count files
    # Args: output path (cmd line input) and name of experiment (cmd line input)
    # Return: directory path

    dir_name = os.path.join(output, exp_name)
    os.system("mkdir -p {}".format(dir_name))

    return dir_name

def filepathList(query, file_type):
    # create filepath from COUNT_LTS, runNumber and fastqFileName
    # Args: queryDB output, the file_type of file to create (default is read_count.tsv and novoalign.log)
    # Return: list of count filepaths

    run_num = list(query['runNumber'])

    base_path = [os.path.join(COUNT_LTS, 'run_{}'.format(x)) for x in run_num]

    fastq_file_basename = [os.path.basename(fileBaseName(x)) for x in query['fastqFileName']]
    read_count_basename = [x + file_type for x in fastq_file_basename]

    if not len(read_count_basename) == len(base_path):
        print('the number of runNumbers and fastqFilePaths is not equal. Please check the query and filter out any lines \
                without run numbers or fastqFileNames as these do not have associated count files')
        sys.exit(1)

    filepath_list = [os.path.join(x, y) for x, y in zip(base_path, read_count_basename)]

    return filepath_list

def moveFiles(file_list, dest_dir, query_len):
    # extract run number and index, cp the files to the destination dir and log the move
    # Args: list of files to be moved, the destination, and the number of rows in the query
    # Return: None

    count = 0
    for file in file_list:
        # throw error/exit if file isn't found in src_dir
        if not os.path.isfile(file):
            print('{} cannot be found and therefore can not be moved. Please check {} for the directory with the \
                   run number in the filename'.format(
                    file, COUNT_LTS))
            sys.exit(1)

        dest_full_path = os.path.join(dest_dir, os.path.basename(file))
        print('...copying {} to {}'.format(os.path.basename(file), dest_full_path))
        cp(file, dest_full_path)
        count = count + 1

    if not count == 2*query_len:
        print("\nThe number of files moved is {}. The number of rows in the query is {}.\n \
               If moving count and bam files (default), the number of files should be twice the number of rows.\n \
               If it is not, Check the query, and {}, and try again".format(count, query_len, COUNT_LTS))
    else:
        print("Your data has been copied to your output directory!")

if __name__ == '__main__':
    main(sys.argv)
