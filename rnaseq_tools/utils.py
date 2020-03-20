import os
import re
import numpy as np
from itertools import combinations, product
import yaml
import sys
import subprocess
import configparser

def decomposeStatus2Bit(status):
    """
    Decompose a value in the column of the output of quality_assess_2 to bit.
    eg 18 = 2 + 16 or the powers of 2 [1.0,4.0]. See templates/qc_config.yaml
    :param status: an int from the status column of the .csv output of quality_assess_2
    :returns: a list of powers of 2 representing the bit
    """
    if status == 0:
        return None
    decomp = []
    for i in np.arange(np.floor(np.log2(status)), -1, -1):
        if status - 2 ** i >= 0:
            decomp.append(i)
            status -= 2 ** i
    return decomp


def makeCombinations(lst):
    """
    Make all possible replicate combinations of the inputted list
    :param lst: a list of items
    :returns: a list of combinations of the inputted list
    """
    if len(lst) < 2:
        return [list]
    combos = []
    for i in range(len(lst), 1, -1):
        for s in combinations(lst, i):
            combos.append(sorted(s))
    return combos


def makeListProduct(lst):
    """
    Make all possible value combination, with each value coming from a each list. Support up to 10 lists.
    :param lst: a list of items
    :returns: all value combinations of the inputted list
    """
    if len(lst) == 0:
        return None
    elif len(lst) == 1:
        return list(product(lst[0]))
    elif len(lst) == 2:
        return list(product(lst[0], lst[1]))
    elif len(lst) == 3:
        return list(product(lst[0], lst[1], lst[2]))
    elif len(lst) == 4:
        return list(product(lst[0], lst[1], lst[2], lst[3]))
    elif len(lst) == 5:
        return list(product(lst[0], lst[1], lst[2], lst[3], lst[4]))
    else:
        sys.exit('ERROR: The number of contrasts is greater than 10. This length beyond what I can handle')


def loadConfig(json_file):
    """
    Load configuration file (JSON) for QC thresholding and scoring
    :param json_file: the config for the rnaseq_pipeline, in particular, is found in templates/qc_config.yaml
    :returns: the read-in json object
    """
    with open(json_file) as json_data:
        d = yaml.safe_load(json_data)
    return d


def mkdirp(path_to_directory):
    """
    Function to create a directory. Equivalent to mkdir -p in bash. Will create directory if the path does not already exist.
    :param path_to_directory: path to a directory to be created if DNE already
    """
    if not os.path.exists(path_to_directory):
        try:
            os.makedirs(path_to_directory)
        except OSError:
            print("Directory {} cannot be created. See mkdir_p function call "
                  "in the script from which this error occurred".format(path_to_directory))

def addForwardSlash(path):
    """
    add forward slash to path
    :param path: a path, typically to the lowest subdir prior to adding a file
    :returns: the path with a trailing slash added
    """
    if not path.endswith('/'):
        path = path + '/'
    return path

def checkCSV(file):
    """
    test whether a given file is a .csv or something else
    :param file: a file that ends with .csv, .xlsx, .tsv, etc.
    :returns: True if .csv, false otherwise
    """
    # test whether a given file is a .csv or .xlsx
    if re.search('\.csv', file):
        return True
    else:
        return False

def fileBaseName(file_name):
    """
    strip everything after first . in path (may be too extreme -- make sure there are no 'necessary' . in filepaths. However, please try to avoid putting a . in a filepath other than the extension(s)
    :param file_name: the base filename of a file (not a path!)
    :returns: the filename stripped of all extensions
    """
    # https://stackoverflow.com/a/46811091
    if '.' in file_name:
        separator_index = file_name.index('.')
        base_name = file_name[:separator_index]
        return base_name
    else:
        return file_name

def pathBaseName(path):
    """
    strips directory structure and extensions. eg /path/to/file.fastq.gz --> file
    :param path: a filepath that includes a directory structure and file (see description)
    :returns: the basename of the file stripped of the directory structure and extensions
    """
    # This gets the basename of a given path, and then strips all file extensions (even if multiple). see fileBaseName
    file_name = os.path.basename(path)
    return fileBaseName(file_name)

def dirName(path):
    """
    get the directory name of the directory one level up from the end of the path provided. eg /path/to/file.fastq.gz returns 'to'
    :param path: a path to a directory or file
    :returns: the name of the directory one level up from the final level of path. see description for example
    """
    if os.path.split(path)[1] == "":
        dirname = os.path.dirname(path)
    else:
        dirname = os.path.basename(path)

    return dirname

def setAttributes(object, expected_attributes, input_dict):
    """
    a makeshift method to mimic a constructor. Check input against list _attributes and notify user if any are not expected
    (a check for entering 'query' instead of 'query_sheet_path'
    :param object: a instance of a class in which you wish to set attributes
    :param expected_attributes: a list of attributes that you expect to have in the class. see StandardData for an example
    :param input_dict: kwargs entered upon instantiation
    """
    for key, value in input_dict.items():
        if key not in expected_attributes:
            print("{} not in expected attributes. This is not a problem unless you expect the data corresponding to "
                      "the values in {} to be automatically standardized or manipulated".format(key, key))
        setattr(object, key, value)

def countsPerMillion(raw_count_path, output_FULL_path):
    """ TODO: re-write this with python subprocess to control input/output of R script
    submit raw_count_path to log2_cpm.R (in tools/)
    :param raw_count_path: path to output of raw_counts.py
    :param output_FULL_path: the full path (including the file and extension) of the output of log2_cpm.R. eg <experiment_name>_log2_cpm.csv
    """
    cmd = 'log2_cpm.R -r {} -o {}'.format(raw_count_path, output_FULL_path)
    executeSubProcess(cmd)

def executeSubProcess(cmd):
    """ TODO: re-do this with package subprocess and store rather than write. write to logger
    executes command, sys.exit with message if the subprocess fails
    :param cmd: the full command to be run
    """
    exit_status = subprocess.call(cmd, shell=True)
    if exit_status == 1:
        sys.exit("{} failed to execute. check the code.".format(cmd))


def configure(object, config_file, config_header):
    """
    reads and sets the attributes in a config_file.ini in type output by configparser.
    NOTE: the object MUST have a attribute object.self_type = 'ObjectName'
    :param object: an object to be configured
    :param config_file: a .ini
    :param config_header: the [header] in the .ini file to read (config format created by configparser)
    """
    # read config file
    config = configparser.ConfigParser()
    config.read(config_file)
    # set attributes for StandardData
    for key, value in config[config_header].items():
        setattr(object, key, value) # by default, values are read in as strings. Currently, all filepaths, so this is good