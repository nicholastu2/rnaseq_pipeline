#!/usr/bin/env python

# Purpose: copy count files from /lts into /scratch/mblab/$USER for analysis
# Input: QueryDB result, location of stored count files, name of experiment, output location (default .)

import os
import sys
import pandas as pd
import re
from move_alignment_count_files import cp
import argparse

COUNT_LTS = '/lts/mblab/Crypto/rnaseq_data/align_expr'

def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query', required=True,
					help='query the database using queryDB.py for the libraries you wish to analyze')
    parser.add_argument('-o', '--output_location', required=True,
					help='Suggested usage: /scratch/$USER or /scratch/$USER/rnaseq  The location where you want to put the directory containing the counts for your analysis')
    parser.add_argument('-n', '--expirement_name', required=True,
						help='The name of this experiment. This will be used as the subdirectory of output')
    return parser.parse_args(argv[1:])

# create experiment directory in output
def createExperimentDir(output, exp_name):
    # create experiment directory to store copied count files
    # Args: output path (cmd line input) and name of experiment (cmd line input)
    # Return: directory path

    dir_name = sys.path.join(output, exp_name)
    os.system("mkdir -p {}".format(dir_name))

    return dir_name

def fastqBasename(file):
    # create list of samples from the fastqFileName column of sample_summary.csv
    # Args: sample_summary.csv (see queryDB in rnaseq_pipeline/tools)
    # Returns: a list of the samples (the fastq file names minus the file extentions)

    # extract column fastqFileName of sample_summary.csv as list
    sample_fastq_list = pd.read_csv(file)['fastqFileName']

    # eliminate file extension (either .fastq.gz or .fastq)
    sample_basename = [re.sub(".fastq.gz|.fastq", "", basename) for basename in sample_fastq_list]

    return sample_basename

# get count filepath
def countFilepathList(query):
    # create filepath from COUNT_LTS, runNumber and fastqFileName
    # Args: queryDB output
    # Return: list of count filepaths

    run_num = list(query['runNumber'])

    base_path = [os.path.join(COUNT_LTS, 'run_{}'.format(x)) for x in run_num]

    fastqFileBasename = [fastqBasename(x) for x in query['fastqFileName']]

    if not len(fastqFileBasename) == len(base_path):
        print('the number of runNumbers and fastqFilePaths is not equal. Please check the query and filter out any lines without run numbers or fastqFileNames as these do not have associated count files')
        sys.exit(1)

    filepath_list = [os.path.join(x,y) for x, y in zip(base_path, fastqFileBasename)]

    return filepath_list

# cp count file to experiment directory
def moveFiles(file_list, dest_dir, query_len):
    # extract run number and index, cp the files to the destination dir and log the move
    # Args: list of files to be moved, the destination, and the number of rows in the query
    # Return: None

    count = 0

    for file in file_list:
            # throw error/exit if file isn't found in src_dir
            if not os.path.isfile(file):
                print('{} cannot be found and therefore can not be moved. Please check {} for the directory with the run number in the filename'.format(file, COUNT_LTS))
                sys.exit(1)

            print('...copying {} to {}'.format(os.path.basename(file), dest_dir))
            cp(file, dest_dir)
        count = count+1

    if not count == query_len:
        print("The number of files moved is {} and the number of rows in the query is {}. Check the query, {}, and try again".format(count, query_len, COUNT_LTS))