#!/usr/bin/env python
import sys
import os
import argparse
import pandas as pd
from rnaseq_tools.CryptoQualAssessAuditObject import CryptoQualAssessAuditObject
from rnaseq_tools.S288C_R54QualAssessAuditObject import S288C_R54QualAssessAuditObject
from rnaseq_tools import utils


# TODO: CURRENTLY ONLY SET UP FOR CRYPTO. NEED TO WRITE S288C_R64QualityAssessmentObject
# TODO: CURRENTLY, ONLY SET UP TO ACCEPT -CC
def main(argv):
    args = parseArgs(argv)
    # parse cmd line arguments and error check paths/values
    print('...parsing cmd line input')
    try:
        if not os.path.isdir(args.align_count_dir):
            raise NotADirectoryError('OutputDirDoesNotExist')
    except NotADirectoryError:
        print('%s does not lead to a valid directory. Check the path and resubmit with working -r' % args.align_count_dir)
    else:
        align_count_path = args.align_count_dir
        output_directory = args.align_count_dir
    try:
        if not os.path.isfile(args.query_sheet_path):
            raise FileNotFoundError('QuerySheetDoesNotExist')
    except FileNotFoundError:
        print('%s does not lead to a valid file. Check and resubmit correct -qs' % args.query_sheet_path)
    except TypeError:
        pass
    else:
        query_sheet_path = args.query_sheet_path

    # get run number if exists for output naming. if DNE, ask user to provide name to insert after run_<using run_num>_summary.csv
    try:
        run_number = utils.getRunNumber(align_count_path)
        # create name for qual_assess
        filename_prefix = 'run_%s' % run_number
    except AttributeError:  # TODO: this will cause a problem if running via batchscript
        filename_prefix = input(
            'No run number detected in input directory name. Enter something to insert in the output directory\n'
            'name: <your_input>_quality_summary.csv: ')
    # store interactive flag
    try:
        interactive_flag = args.interactive
    except AttributeError:
        interactive_flag = False

    # read in query sheet # TODO: GENERALIZE THIS INTO EITHER STANDARDDATA OR UTILS. RETURN AS DICT. DO THIS AFTER ADDING ORGANISM COLUMN TO METADATA SPECS
    query_df = utils.readInDataframe(query_sheet_path)
    query_fastq_list = list(query_df.fastqFileName)

    # extract bam file names
    bam_list = utils.extractFiles(align_count_path, '.bam')
    # filter bam_list for files in the query sheet
    filtered_bam_list = [x for x in bam_list if os.path.basename(x).replace('_sorted_aligned_reads_with_annote.bam', '.fastq.gz') in query_fastq_list]
    # extract novoalign logs
    novoalign_logs = utils.extractFiles(align_count_path, 'novoalign.log')
    filtered_novoalign_logs = [x for x in novoalign_logs if os.path.basename(x).replace('_novoalign.log', '.fastq.gz') in query_fastq_list]
    # extract count file list
    count_list = utils.extractFiles(align_count_path, 'read_count.tsv')
    filtered_count_list =  [x for x in count_list if os.path.basename(x).replace('_read_count.tsv', '.fastq.gz') in query_fastq_list]
    # from count_list, get convert to a list of fastq.gz names
    extracted_sample_fastq_list = [os.path.basename(x.replace('_read_count.tsv', '.fastq.gz')) for x in count_list]
    if len(filtered_bam_list) != len(filtered_count_list) or len(filtered_bam_list) != len(filtered_novoalign_logs):
        sys.exit('The number of bam_files, count_files and/or log_files does not match. Check file contents')

    # all crypto records will have genotype beginning with CNAG_
    crypto_query_df = query_df[~query_df.genotype1.isna() & query_df.genotype1.str.startswith('CNAG') & query_df.fastqFileName.isin(extracted_sample_fastq_list)]
    yeast_query_df = query_df[(~(query_df.genotype1.isna() | query_df.fastqFileName.isin(crypto_query_df.fastqFileName)) & query_df.fastqFileName.isin(extracted_sample_fastq_list))]

    # create list to store qual_assess dataframes
    qual_assess_df_list = []

    if len(crypto_query_df) > 0:
        # if coverage_check is passed in cmd line, include query and coverage_check_flag in constructor (automatically sets some values #TODO make this a function with arugmnets to pass so as not to repeat entire constructor)
        print('...compiling KN99 samples information')
        crypto_qa_object = CryptoQualAssessAuditObject(organism = 'KN99',
                                                       bam_file_list=filtered_bam_list,
                                                       count_file_list=filtered_count_list,
                                                       novoalign_log_list=filtered_novoalign_logs,
                                                       coverage_check_flag=True,
                                                       query_df=crypto_query_df,
                                                       config_file=args.config_file,
                                                       interactive=interactive_flag)

        # add dataframe to list
        try:
            qual_assess_df_list.append(crypto_qa_object.qual_assess_df)
        except AttributeError:
            error_msg = 'There was an error appending the KN99 qual assess dataframe. Check the paths in the query sheet and align_counts directory'
            crypto_qa_object.logger.debug(error_msg)
            print(error_msg)

    if len(yeast_query_df) > 0:
        yeast_qa_object = S288C_R54QualAssessAuditObject(organism='S288C_R64',
                                                           bam_file_list=filtered_bam_list,
                                                           count_file_list=filtered_count_list,
                                                           novoalign_log_list=filtered_novoalign_logs,
                                                           query_path=args.query_sheet_path,
                                                           config_file=args.config_file,
                                                           interactive=interactive_flag)
        print('...compiling S288C_R64 alignment information')
        # create dataframes storing the relevant alignment and count metadata from the novoalign and htseq logs
        try:
            qual_assess_df_list.append(yeast_qa_object.qual_assess_df)
        except AttributeError:
            error_msg = 'There was an error appending the S288C_R64 qual assess dataframe. Check the paths in the query sheet and align_counts directory'
            crypto_qa_object.logger.debug(error_msg)
            print(error_msg)

    # combine dataframes, if both organisms present
    print('...creating quality_assessment sheet for %s' %filename_prefix)
    combined_qual_assess_1_df = pd.concat(qual_assess_df_list)

    # create filename
    quality_assessment_filename = "%s_sequence_quality_summary.csv" % filename_prefix
    output_path = os.path.join(output_directory, quality_assessment_filename)
    print('writing output to %s' % output_path)
    combined_qual_assess_1_df.to_csv(output_path, index=False)


def parseArgs(argv):
    parser = argparse.ArgumentParser(description="This script summarizes the output from pipeline wrapper.")
    parser.add_argument("-ac", "--align_count_dir", required=True,
                        help="[REQUIRED] Directory with files in the following subdirectories: align, count, logs. Output from raw_count.py and log2cpm.R must be in count directory.")
    parser.add_argument("-qs", "--query_sheet_path",
                        help="[REQUIRED] Path to query sheet")
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
