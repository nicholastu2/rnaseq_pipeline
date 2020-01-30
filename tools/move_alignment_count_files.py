mport glob
import re
import os
import sys
import argparse
from shutil import copy2 as cp


def main(argv):
    # read in cmd line args
    args = parseArgs(argv)

    alignment_files = getAlignmentFiles(args.alignment_directory)
    log_files = getLogFiles(args.log_directory)

    moveFiles(alignment_files, args.destination_path, args.log)
    moveFiles(log_files, args.destination_path, args.log)

    count_sheet = countSheet(args.log)

    cp(count_sheet, args.metadata_count_dir)

def parseArgs(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-ad', '--alignment_directory', required = True,
                        help = 'path .bam and .tsv files (typically alignment/<some sub dir>)')
    parser.add_argument('-ld', '--log_directory', required = True,
                        help = 'path to htseq, novoalign and novosort logs')
    parser.add_argument('-dp', '--destination_path', required = True,
                        help = 'path to destination directory')
    parser.add_argument('-l', '--log', required = True,
                        help = 'full path (path and file name) to a .tsv log file. eg) logs/move_4010_samples.tsv')

    return parser.parse_args(argv[1:])

def addForwardSlash(path):
    # adds a forward slash to a directory path if one is not already present
    # Args: a path
    # Returns: the path with a forward slash (does not add one if already there)

    if not path[-1] == '/':
        path = path + '/'
    return path

def getFileNames(dir, *args):
    # return file names with certain suffix
    # Args: the directory to search (dir) and a list of suffixes to search for in the dir
    # Return: file names as a list
    all_files = []

    for suffix in args:
        dir = addForwardSlash(dir)
        search_pattern = dir + '*{}'.format(suffix)
        filepath_list = glob.glob(search_pattern)
        all_files.append(filepath_list)

    return filepath_list

def getAlignmentFiles(alignment_dir):
    # get list of alignment files
    # args alignment_dir path
    # Return: list of files with suffixes listed below

    suffixes = ['_aligned_reads.bam', '_read_count.tsv']

    alignment_files = getFileNames(alignment_dir, suffixes)

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
    # Return: None. log is written where indicated

    with open(log_file, 'a+') as cp_log:
        for file in file_list:

            # TODO: point of weakness -- IMPROVE THE REGEX/RUN_NUM EXTRACTION
            regex_num = r"(\d*)_.*.*"
            run_num = re.search(regex_num, file).group(1)
            # TODO: SEE ABOVE
            regex_index = r".*_Index\d*_([ATGC]*)"
            index = re.search(regex_index, file).group(1)

            destination_dir = addForwardSlash(destination_dir)
            destination_file_path = os.path.join(destination_dir, os.path.basename(file))

            if not os.path.isfile(file):
                print('{} cannot be found and therefore can not be moved. Please check cmd line inputs')
                sys.exit(1)
            else:
                cp(file, destination_file_path)
                cp_log.write('{}\t{}\t{}\t{}\n'.format(run_num, index, file, destination_file_path))
                print('...moving {} to {}'.format(os.path.basename(file), destination_dir))

def countSheet(log_file):
    # parses log_file for read_counts.tsv
    # Args: the completed log_file (must be after moves)
    # Returns: the filepath to the new count sheet



if __name__ == '__main__':
    main(sys.argv)

