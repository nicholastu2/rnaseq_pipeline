#!/usr/bin/env python
import sys
import os
import argparse
from rnaseq_tools.CryptoQualityAssessmentObject import CryptoQualityAssessmentObject
from rnaseq_tools.S288C_R64QualityAssessmentObject import S288C_R64QualityAssessmentObject
from rnaseq_tools import utils


# TODO: CURRENTLY ONLY SET UP FOR CRYPTO. NEED TO WRITE S288C_R64QualityAssessmentObject
# TODO: CURRENTLY, ONLY SET UP TO ACCEPT -CC
def main(argv):
    args = parseArgs(argv)
    # parse cmd line arguments and error check paths/values
    print('...parsing cmd line input')
    try:
        if not args.organism in ['KN99', 'S288C_R64', 'H99']:
            raise ValueError('InvalidOrganism: %s' %args.organism)
        organism=args.organism
    except ValueError:
        print('The organism must be KN99, S288C_R64, or H99 (use that spelling/capitalization exactly)')
    try:
        if not os.path.isdir(args.align_count_dir):
            raise NotADirectoryError('OutputDirDoesNotExist')
    except NotADirectoryError:
        print('%s does not lead to a valid directory. Check the path and resubmit with working -r' % args.align_count_dir)
    else:
        align_count_path = args.align_count_dir
    try:
        if not os.path.isdir(args.output_dir):
            raise NotADirectoryError
    except NotADirectoryError:
        print('%s does not lead to a valid directory. check the path and resubmit with correct -o' % args.output_dir)
    else:
        output_directory = args.output_dir
    try:
        if not os.path.isfile(args.query_sheet_path):
            raise FileNotFoundError('QuerySheetDoesNotExist')
    except FileNotFoundError:
        print('%s does not lead to a valid file. Check and resubmit correct -qs' % args.query_sheet_path)
    except TypeError:
        pass

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

    # extract bam file names
    bam_list = utils.extractFiles(align_count_path, '.bam')
    # extract novoalign logs
    novoalign_logs = utils.extractFiles(align_count_path, 'novoalign.log')
    # extract count file list
    count_list = utils.extractFiles(align_count_path, 'read_count.tsv')
    if len(bam_list) != len(count_list) or len(bam_list) != len(novoalign_logs):
        sys.exit('The number of bam_files, count_files and/or log_files does not match. Check file contents')

    # create filename
    quality_assessment_filename = "%s_quality_summary.csv" % filename_prefix
    output_path = os.path.join(output_directory, quality_assessment_filename)

    if organism == 'KN99':
        # if coverage_check is passed in cmd line, include query and coverage_check_flag in constructor (automatically sets some values #TODO make this a function with arugmnets to pass so as not to repeat entire constructor)
        crypto_qa = CryptoQualityAssessmentObject(bam_file_list=bam_list,
                                                  count_file_list=count_list,
                                                  novoalign_log_list=novoalign_logs,
                                                  coverage_check_flag=True,
                                                  query_path=args.query_sheet_path,
                                                  config_file=args.config_file,
                                                  interactive=interactive_flag)

        print('...compiling alignment information')
        # create dataframes storing the relevant alignment and count metadata from the novoalign and htseq logs
        crypto_qual_assess_1_df = crypto_qa.compileAlignCountMetadata()

        print('writing output to %s' % output_path)
        crypto_qual_assess_1_df.to_csv(output_path, index=False)

    elif organism == 'S288C_R64':
        yeast_qa = S288C_R64QualityAssessmentObject(bam_file_list=bam_list,
                                                    count_file_list=count_list,
                                                    novoalign_log_list=novoalign_logs,
                                                    query_path=args.query_sheet_path,
                                                    config_file=args.config_file,
                                                    interactive=interactive_flag)
        print('...compiling alignment information')
        # create dataframes storing the relevant alignment and count metadata from the novoalign and htseq logs
        s288c_r64_qual_assess_1_df = yeast_qa.compileAlignCountMetadata()

        print('writing output to %s' % output_path)
        s288c_r64_qual_assess_1_df.to_csv(output_path, index=False)


def parseArgs(argv):
    parser = argparse.ArgumentParser(description="This script summarizes the output from pipeline wrapper.")
    parser.add_argument("-ac", "--align_count_dir", required=True,
                        help="[REQUIRED] Directory for alignment log (novoalign) and count (htseq) files.")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="[REQUIRED] File path to the directory you wish to deposit the summary. Note: the summary will be called run_###_summary.csv")
    parser.add_argument("-qs", "--query_sheet_path",
                        help="[REQUIRED] Path to query sheet filtered for the files contained in the path passed to -r")
    parser.add_argument("-g", "--organism",
                        help="[REQUIRED] Either KN99, S288C_R64, or H99")
    parser.add_argument("-pc", "--perturbation_check", action='store_true',
                        help="[OPTIONAL] For Crypto experiments. Set this flag to add coverage and overexpression columns. Note: this makes the script take a long time to complete")
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
