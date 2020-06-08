#!/usr/bin/env python3
"""
   For a given alignment/count job, create a nextflow config file. See templates for example.
   Currently, the expectation is that you run this on the cluster.
   create_nextflow_config.py -qs /path/to/database/sheet/describing/runs/to/process
"""
import sys
import argparse
import os
import numpy as np
import pandas as pd
from rnaseq_tools.DatabaseObject import DatabaseObject
from rnaseq_tools import utils


def main(argv):
    """ main method
    :param argv: cmd line arguments
    """
    # parse cmd line arguments
    args = parseArgs(argv)
    query_sheet_path = args.query_sheet
    try:
        if not os.path.isfile(query_sheet_path):
            raise FileNotFoundError
    except FileNotFoundError:
        print('Query sheet path not valid. Check and try again.')
    try:
        interactive_flag = args.interactive
    except AttributeError:
        interactive_flag = False

    # instantiate DatabaseObject --> mostly this will be for access to StandardData paths
    db = DatabaseObject(query_sheet_path=query_sheet_path, config_file=args.config_file, interactive=interactive_flag)
    # read in dataframe
    db.query_df = utils.readInDataframe(db.query_sheet_path)
    db.query_df['organism'] = np.where(db.query_df['genotype'].str.startswith('CNAG'), 'KN99', 'S288C_R64')
    db.query_df['libraryDate'] = pd.to_datetime(db.query_df['libraryDate'])
    db.query_df['strandedness'] = np.where(db.query_df['libraryDate'] > '2015-10-25', 'reverse', 'no')
    # add leading zero to runNumber, if necessary
    db.query_df['runNumber'] = db.query_df['runNumber'].astype(str)
    # new dictionary to store run_directory in dataframe
    run_directory_list = []
    for index, row in db.query_df.iterrows():
        if int(row['runNumber']) in (db._run_numbers_with_zeros):
            run_number = str(db._run_numbers_with_zeros[int(row['runNumber'])])
        else:
            run_number = str(row['runNumber'])
        run_directory = 'run_' + run_number + '_samples'
        run_directory_list.append(run_directory)
        fastq_filename = os.path.basename(row['fastqFileName'])
        fastq_fullpath = os.path.join(db.lts_sequence, run_directory, fastq_filename)
        db.query_df.loc[index, 'fastqFileName'] = fastq_fullpath
    db.query_df['runDirectory'] = run_directory_list

    # filter
    nextflow_fastqfile_df = db.query_df[['runDirectory', 'fastqFileName', 'organism', 'strandedness']]
    print(nextflow_fastqfile_df)
    # write out
    output_path = os.path.join(db.job_scripts, 'nextflow_fastqfile_list'+ '_' + db.year_month_day + '_' + utils.hourMinuteSecond()+'.csv')
    print('writing out to %s' %output_path)
    nextflow_fastqfile_df.to_csv(output_path, index=False)



def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="For a given alignment/count job, create a nextflow config file")
    parser.add_argument("-qs", "--query_sheet", required=True,
                        help='[REQUIRED] A .csv subset of metadata database describing the set of files you wish to process.\n'
                             'These will be grouped by run number')
    parser.add_argument('--config_file', default='/see/standard/data/invalid/filepath/set/to/default',
                        help="[OPTIONAL] default is already configured to handle the invalid default path above in StandardDataObject.\n"
                             "Use this flag to replace that config file")
    parser.add_argument('--interactive', action='store_true',
                        help="[OPTIONAL] set this flag (only --interactive, no input necessary) to tell StandardDataObject not\n"
                             "to attempt to look in /lts if on a compute node on the cluster")

    args = parser.parse_args(argv[1:])
    return args

if __name__ == "__main__":
    main(sys.argv)