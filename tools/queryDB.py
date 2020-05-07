#!/usr/bin/env python

# name: queryDB
# purpose: parse rnaseq metadata
# input: bio sample database main directory; filepath for script output; json with search terms; name of query
# output: four .csv : expr.lookup.txt has paths to the count matricies;
#                     fastq.lookup.txt has paths to the fastqs;
#                     queriedDB.csv is the complete database joined from input directory;
#                     sample_summary.csv is the filtered table for input into rnaseq_pipe
# written by: implemented in python by chase mateusiak(chase.mateusiak@gmail.com), based on sanji bhavsar code in R.
#             See https://github.com/sanjibhavsar/rnaseq_database for original work
# date included in rnaseq_pipe: 1/20/2020
# See the end of the script for a description of the environment used to write/test

import os
import sys
import argparse
import time

# list of subdirectories of datadir to be used as keys of dictionary. When these files are searched, only .csv and .xlsx
# are listed for concatenating. The search through these directories is not recursive.
# this is not a user input b/c there may be subdirectories that we do not wish to search.
# HOWEVER, the functions that require datadir_keys take it as input, which allows those functions to be used to search
# different subdirs as needed
datadir_keys = ['fastqFiles', 'library', 's2cDNASample', 's1cDNASample', 'rnaSample', 'bioSample']

def main(argv):
    # read in cmd line args
    args = parseArgs(argv)

    # get filepaths of all sheets in the various subdirs of the datadir you passed in cmd line
    datadir_dict = getFilePaths(args.database)
    # combine on the common columns the files in the subdirs (the datadir_keys) passed in cmd line
    combined_df = createDB(datadir_dict)

    # if there is a query, parse and write query
    if args.json:
        # query the combined db based on json input from cmd line
        query_df = queryDB(combined_df, args.json)
        # write out
        query_name = os.path.basename(args.json).split('.')[0]
        query_output = os.path.join(args.output, query_name)
        query_df.to_csv(query_output + '.csv', index=False)

    # if user enters -pf, print full database
    if args.print_full:
        time_stamp = (time.strftime("%Y_%m_%d", time.gmtime()))
        combined_output = os.path.join(args.output, 'combined_df_{}.csv'.format(time_stamp))
        combined_df.to_csv(combined_output, index=False)

def parseArgs(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--database', required = True,
                        help = 'topmost directory of metadata database. On cluster, /scratch/mblab/database-files. \
                                If using the rnaseq_pipeline module, you may use $METADATA. Do make sure that $METADATA \
                                is up to date by running git pull in the directory')
    parser.add_argument('-j', '--json',
                        help = 'path to json file used to parse metadata. See ')
    parser.add_argument('-o', '--output', required = True,
                        help = 'filepath to directory to intended queryDB output')
    parser.add_argument( '-pf', '--print_full', action='store_true',
                         help = 'Use this in the absence of -j to print out the full metadata database. The name will be combined_df_[date].csv. \
                         Use it in addition to -j to print out both the query and the full database. Note: simply add -pf. No value is necessary')

    return parser.parse_args(argv[1:])



if __name__ == '__main__':
	main(sys.argv)
