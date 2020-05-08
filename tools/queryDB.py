#!/usr/bin/env python
"""
    A script to operate DatabaseObject in rnaseq_tools. If you need more flexibility, do the following in htcf
        cd /scratch/mblab/<user>
        ml rnaseq_pipeline
        python
        #>>> from rnaseq_tools.DatabaseObject import DatabaseObject
        #>>> database_object = DatabaseObject() <-- see rnaseq_tools.DatabaseObject script for help on usage
    written by: implemented in python by chase mateusiak(chase.mateusiak@gmail.com), based on sanji bhavsar code in R.
                 See https://github.com/sanjibhavsar/rnaseq_database for original work
    date included in rnaseq_pipe: 1/20/2020
    revised: 5/7/2020

    from cluster (REMEMBER that this WILL create a directory rnaseq_pipeline in /scratch/mblab/<user>)
    usage: queryDB.py -pf # this will print the full database from /scratch/mblab/<user>/database_files
"""

import os
import sys
import argparse
from rnaseq_tools import utils
from rnaseq_tools.DatabaseObject import DatabaseObject
from rnaseq_tools.StandardDataObject import StandardData

def main(argv):
    # read in cmd line args
    args = parseArgs(argv)
    print('...parsing arguments')
    # if not cluster is passed, create logger in PWD
    if args.not_cluster:
        logger_path = os.path.join(os.getcwd(), 'queryDB_log_%s.log' %utils.yearMonthDay())
        logger = utils.createLogger(logger_path, 'queryDB.py')
    # if -nc is not passed, assume on HTCF and call StandardData to check/load rnaseq_pipeline file structure and data
    else:
        sd = StandardData()
        sd.standardDirectoryStructure()
        sd.setStandardDataLogger()
        logger_path = sd.log_file
        logger = utils.createLogger(logger_path, __name__)
    print('queryDB log can be found at: %s' %logger_path)
    logger.debug('cmd line arguments are: %s'% args)

    # read in and check cmd line arguments
    database_path = args.database
    if database_path is not None and not os.path.exists(database_path):
        raise FileNotFoundError('DatabaseFileDoesNotExist')
    elif database_path is None:
        database_path = sd.scratch_database_files
    filter_json_path = args.json
    if filter_json_path is not None and not os.path.isfile(filter_json_path):
        raise FileNotFoundError('QueryJsonDoesNotExist')
    output_directory = args.output_directory
    if output_directory is not None and not os.path.exists(output_directory):
        raise FileNotFoundError('OutputDirectoryDoesNotExist')
    print('...compiling database')
    # create database object (see rnaseq_tools.DatabaseObject for description of all functions
    database_object = DatabaseObject(database_path, logger_path=logger_path, filter_json_path=filter_json_path)
    database_object.setDatabaseDataframe()

    # filter database and print to output_directory, if json is present
    if database_object.filter_json_path is not None:
        print('...filtering database')
        database_object.filterDatabaseDataframe()
        output_filename = utils.pathBaseName(database_object.filter_json_path)
        filtered_output_path = os.path.join(output_directory, output_filename)
        print('printing filtered database to: %s' %filtered_output_path)
        database_object.filtered_database_df.to_csv(filtered_output_path, index=False)

    # if user enters -pf, print full database
    if args.print_full:
        year_month_day = utils.yearMonthDay()
        full_database_output_path = os.path.join(output_directory, 'combined_df_{}.csv'.format(year_month_day))
        print('printing full database to: %s' %full_database_output_path)
        database_object.database_df.to_csv(full_database_output_path, index=False)

def parseArgs(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_directory', required=True,
                        help='[REQUIRED] Filepath to directory to intended queryDB output')
    parser.add_argument('-nc', '--not_cluster', action='store_true',
                        help='[OPTIONAL] Set flag (no input) to use off HTCF')
    parser.add_argument('-d', '--database', default=None,
                        help='[OPTIONAL] Default is database_files in /scratch/mblab/user/rnaseq_pipeline. '
                             'If entered, use topmost directory of metadata database. On cluster, /scratch/mblab/database-files.')
    parser.add_argument('-j', '--json', default=None,
                        help='[OPTIONAL] Path to json file used to parse metadata. See ')
    parser.add_argument('-pf', '--print_full', action='store_true',
                        help='[OPTIONAL] Use this in the absence of -j to print out the full metadata database. The name will be combined_df_[date].csv. \
                         Use it in addition to -j to print out both the query and the full database. Note: simply add -pf. No value is necessary')

    return parser.parse_args(argv[1:])


if __name__ == '__main__':
    main(sys.argv)
