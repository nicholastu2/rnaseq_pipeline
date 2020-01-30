#!/usr/bin/env python
import glob
import re
import os
import sys
import argparse
import pandas as pd
from shutil import copy2 as cp

def main(argv):
    # read in cmd line args
    args = parseArgs(argv)

    # get file_paths
    alignment_files = getAlignmentFiles(args.alignment_directory)
    log_files = getLogFiles(args.log_directory)

    align_dict = moveFiles(alignment_files, args.destination_path, args.report)
    log_dict = moveFiles(log_files, args.destination_path, args.report)

    writeCountSheet(align_dict, log_dict, args.count_metadata)

def parseArgs(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-ad', '--alignment_directory', required = True,
                        help = 'path .bam and .tsv files (typically alignment/<some sub dir>)')
    parser.add_argument('-ld', '--log_directory', required = True,
                        help = 'path to htseq, novoalign and novosort logs')
    parser.add_argument('-dp', '--destination_path', required = True,
                        help = 'path to destination directory')
    parser.add_argument('-r', '--report', required = True,
                        help = 'full path (path and file name) to a .tsv report file. This is redundant if the move is successful. Store in scratch folder meant for non-essential logs/repors')
    parser.add_argument('-c', '--count_metadata', required = True,
                        help = 'path to alignCount directory in metadata database')

    return parser.parse_args(argv[1:])

def writeCountSheet(align_dict, log_dict, count_metadata):

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

def moveFiles(file_list, destination_dir, log_file):
    # extract run number and index, cp the files to the destination dir and log the move
    # Args: list of files to be moved (must have file paths from PWD, destination diretory and a log file (complete path including name of file)
    # Return: dictionary of key (run_num + index) : values [list destination file names]
    # Actions: Unless the dry-run flag is passed, this funciton copies src files to destination path

    move_dict={}

    with open(log_file, 'a+') as cp_log:
        for file in file_list:
            try:
                regex_num = r"(\d*)_.*.*"
                run_num = re.search(regex_num, file).group(1)
            except AttributeError:
                print('no run_number found in file: {}'.format(file))
            try:
                regex_index = r".*_Index\d*_([ATGC]*)"
                index = re.search(regex_index, file).group(1)
            except AttributeError:
                print('no index found in file {}'.format(file))

            # create destination_file_path
            destination_dir = addForwardSlash(destination_dir)
            destination_file_path = os.path.join(destination_dir, os.path.basename(file))

            # throw error/exit if file isn't found in src_dir
            if not os.path.isfile(file):
                print('{} cannot be found and therefore can not be moved. Please check cmd line inputs')
                sys.exit(1)

            print('...moving {} to {}'.format(os.path.basename(file), destination_dir))
            cp(file, destination_file_path)
            cp_log.write('{}\t{}\t{}\t{}\n'.format(run_num, index, file, destination_file_path))
            move_dict.setdefault((run_num, index), []).append(destination_file_path)

    return move_dict


if __name__ == '__main__':
    main(sys.argv)
