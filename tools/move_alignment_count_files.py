#!/usr/bin/env python
"""
    move alignment_count files from user scratch space to long term storage (or anywhere else for that matter, but this is the intent)

    from the directory structure created by StandardDataObject.standardDirectoryStructure()
    usage: move_alignment_count_files.py -r reports/run_1234 -rn 1234 -d lts_align_expr
"""
import glob
import os
import sys
import argparse
import pandas as pd
import time
from rnaseq_tools import utils

# TODO: use StandardData and incorporate logging

def main(argv):

    # read in cmd line args
    args = parseArgs(argv)

    # get file_paths to alignment and count files
    alignment_files = getAlignmentFiles(args.reports)
    novo_log_files = getLogFiles(args.reports)

    # create destination directory
    dest_path_complete = os.path.join(args.destination_path, 'run_{}'.format(args.run_number))
    utils.mkdirp(dest_path_complete)  # this will not overwrite directory if already exists

    # move the alignment files
    moveFiles(alignment_files, args.destination_path, args.run_number)
    # move the log files
    moveFiles(novo_log_files, args.destination_path, args.run_number)

    # move the pipeline info directory. If a pipeline info directory already exists, append date time to pipeline_info dir in question
    pipeline_info = os.path.join(args.reports, 'pipeline_info')
    if os.path.isdir(pipeline_info):
        if os.path.isdir(os.path.join(dest_path_complete, 'pipeline_info')):
            os.system('rsync -aHv {} {}'.format(pipeline_info, os.path.join(dest_path_complete, "pipeline_info" + '_' + time.strftime("%Y%m%d-%H%M%S"))))
            #shutil.copytree(pipeline_info, os.path.join(dest_path_complete, "pipeline_info" + '_' + time.strftime("%Y%m%d-%H%M%S")))
        else:
            os.system('rsync -aHv {} {}'.format(pipeline_info, os.path.join(dest_path_complete, "pipeline_info")))
            #shutil.copytree(pipeline_info, os.path.join(dest_path_complete, "pipeline_info"))

def parseArgs(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reports', required = True,
                        help='[REQUIRED] the reports/run_## directory in which the alignment and count files are')
    parser.add_argument('-rn', '--run_number', required=True,
                        help = '[REQUIRED] The run number corresponding to the set of fastq files')
    parser.add_argument('-d', '--destination_path', required = True,
                        help = '[REQUIRED] Suggested usage: /lts/mblab/Crypto/rnaseq_data/lts_align_expr   path to destination directory')

    return parser.parse_args(argv[1:])

def writeMoveLog(align_dict, log_dict, count_metadata):
    # This generates the alignCounts automatically generated metadata sheet for the database-files
    # Args: a dictionry of alignment files, log files, and the path to database-files/alignCount
    # Output: a dataframe written to the path provided in the cmd line input -c

    # combine dictionaries
    for key, value in log_dict.items():
        align_dict.setdefault(key, []).append(value)
    # flatten the list entries for each key and rename
    count_dict = {key: flattenList(value) for key, value in align_dict.items()}
    df = pd.DataFrame.from_dict(count_dict, orient='index')
    df[['runNumber', 'index1Sequence']] = pd.DataFrame(df.index.to_list(), index=df.index)
    df = df.rename(columns={0: 'sortedAlignment', 1: "readCounts", 2: "htseqLog", 3: "novoalignLog", 4: "novosortLog"})
    df = df.reset_index()
    df = df.drop(columns='index')
    cols_to_move = ['runNumber', 'index1Sequence']
    df = df[cols_to_move + [col for col in df.columns if col not in cols_to_move]]
    meta_dir = os.path.join(count_metadata,'run_{}.csv'.format(df.at[0,'runNumber']))
    df.to_csv(meta_dir, index=False)

def flattenList(list_lists):
    # flattens a list with list elements
    # Args: a list with elements that have either all list terms or list + not list terms
    # Returns: flattened list

    flat_list = []
    for item in list_lists:
        if isinstance(item, list):
            for sub_item in item:
                flat_list.append(sub_item)
        else:
            flat_list.append(item)
    return flat_list

def addForwardSlash(path):
    # adds a forward slash to a directory path if one is not already present
    # Args: a path
    # Returns: the path with a forward slash (does not add one if already there)

    if not path[-1] == '/':
        path = path + '/'
    return path

def getFileNames(dir, list_args):
    # return file names with certain suffix
    # Args: the directory to search (dir) and a list of suffixes to search for in the dir
    # Return: file names as a list
    all_files = []

    for suffix in list_args:
        dir = addForwardSlash(dir)
        search_pattern = dir + '*{}'.format(suffix)
        filepath_list = glob.glob(search_pattern)
        all_files = all_files + filepath_list

    return all_files

def getAlignmentFiles(alignment_dir):
    # get list of alignment files
    # args alignment_dir path
    # Return: list of files with suffixes listed below

    suffixes = ['_aligned_reads.bam', '_read_count.tsv']

    alignment_files = []
    all_alignment_files = getFileNames(alignment_dir, suffixes)
    # remove non sorted aligned reads
    for file in all_alignment_files:
        if 'aligned_reads' in file:
            if 'sorted' in file:
                alignment_files.append(file)
            else:
                pass
        else:
            alignment_files.append(file)

    return alignment_files

def getLogFiles(log_dir):
    # get list of file names in log dir to be moved
    # Args: log dir path
    # Returns: list of log files to be moved

    suffixes = ['_htseq.log', '_novoalign.log', '_novosort.log']

    htseq_log_files = getFileNames(log_dir, suffixes)

    return htseq_log_files

def moveFiles(file_list, destination_dir, run_num):
    """
        extract run number and index, cp the files to the destination dir and log the move
        :param file_list: list of files
        :param destination_dir: destination directory
        :param run_num: run number of file
        :returns: dictionary of key (run_num + index) : values [list destination file names]
        :actions: Unless the dry-run flag is passed, this funciton copies src files to destination path
    """

    for file in file_list:

        # create destination_file_path. This will be in pattern inputted_dest_path/run_####/.*.bam etc
        run_dir = 'run_{}'.format(run_num)
        destination_dir = addForwardSlash(destination_dir)
        destination_path_intermediate = os.path.join(destination_dir, run_dir)
        destination_file_path = os.path.join(destination_path_intermediate, os.path.basename(file))

        # throw error/exit if file isn't found in src_dir
        if not os.path.isfile(file):
            print('{} cannot be found and therefore can not be moved. Please check cmd line inputs'.format(file))
            sys.exit(1)

        print('...moving {} to {}'.format(os.path.basename(file), destination_path_intermediate))
        cmd = 'rsync -aHv {} {}'.format(file, destination_file_path)
        utils.executeSubProcess(cmd)


if __name__ == '__main__':
    main(sys.argv)
