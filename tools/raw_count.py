#!/usr/bin/env python
import sys
import argparse
import os
import glob
from rnaseq_tools.OrganismDataObject import OrganismData
from rnaseq_tools import utils

# TODO: Update with OrganismData object (no more need to input gene list). Better commeting and explanation of each step
#  in createCountMatrix. Try to abstract this for use when making any matrix (counts qual_assess_1 and 2, etc)
# TODO: for runs/directories with more than one organism, create two sheets (input query_sheet for this)

# NOTE: this is set up very specifically assuming that the counts are output from the nextflow pipeline

def main(argv):

    args = parseArgs(argv)

    try:
        if not os.path.isdir(args.count_directory):
            raise NotADirectoryError('ERROR: %s does not exist.' % args.count_directory)
        if not os.path.isfile(args.query_sheet):
            raise NotADirectoryError('ERROR: %s does not exist.' % args.count_directory)
    except FileNotFoundError:
        print('path to %s does not exist')
    else:
        count_dirpath = args.count_directory
        query_sheet_path = args.query_sheet
        query_df = utils.readInDataframe(query_sheet_path)

    # extract count files from count_dir
    count_dir_file_list = glob.glob(os.path.join(count_dirpath, '*read_count.tsv'))

    # TODO: SOME ERROR CHECKING ON THE FASTQFILENAME?
    # all crypto records will have genotype beginning with CNAG_, used this to extract list of crypto and yeast samples from query
    crypto_sample_list = list(query_df[query_df.genotype1.str.startswith('CNAG')].fastqFileName) #TODO: after metadata organism column added, update this section
    s288c_r64_sample_list = list(query_df[~query_df.genotype1.str.startswith('CNAG')].fastqFileName)

    # split list of count files based on membership in dataframes above
    count_files_by_organism_dict = {'KN99': [x for x in count_dir_file_list if os.path.basename(x.replace('_read_count.tsv', '.fastq.gz')) in crypto_sample_list],
                                    'S288C_R64': [x for x in count_dir_file_list if os.path.basename(x.replace('_read_count.tsv', '.fastq.gz')) in s288c_r64_sample_list]}

    # create and write out count sheets
    for organism, count_file_list in count_files_by_organism_dict.items():
        if len(count_file_list) > 0:
            od = OrganismData(organism = organism, config_file=args.config_file, interactive=args.interactive)
            count_df = od.createCountSheet(count_file_list)
            output_path = os.path.join(utils.dirPath(utils.dirPath(count_file_list[0])), '%s_raw_count.csv' %organism)
            print('writing count file to %s' %output_path)
            count_df.to_csv(output_path, index=False)

def parseArgs(argv):
    parser = argparse.ArgumentParser(description="Create raw counts for each organism in a count directory.\n"
                                                 "To perform this on more than one directory, a batchscript should be\n"
                                                 "used.")
    parser.add_argument('-c', '--count_directory', required=True,
                        help='[REQUIRED] The directory with has subdirectories align, count and logs, created by align/count nextflow pipeline.\n'
                             ' The raw_count csv will be output in this directory with the name run_number/<organism>_raw_count.csv')
    parser.add_argument('-qs', '--query_sheet', required=True,
                        help='[REQUIRED] path to query sheet containing at least the samples in the count directory')
    parser.add_argument('--config_file', default='/see/standard/data/invalid/filepath/set/to/default',
                        help="[OPTIONAL] default is already configured to handle the invalid default path above in StandardDataObject.\n"
                             "Use this flag to replace that config file. Note: this is for StandardData, not nextflow")
    parser.add_argument('--interactive', action='store_true',
                        help="[OPTIONAL] This flag functions to disable StandardData from attempting to soft link to /lts on htcf."
                             " Default is True, which is the setting for running this script in an interactive session. You should run"
                             " this is an interactive session. To do so, simply enter the command interactive and hit enter before running"
                             " this script. To extract genome files or soft link to /lts, set this flag to False. Just make sure you are not in"
                             " an interactive session, which again, you generally should be to run this. The sys admins will thank you.")
    return parser.parse_args(argv[1:])


if __name__ == '__main__':
    main(sys.argv)

