#!/usr/bin/env python
import sys
import os
import argparse
import pandas as pd
from rnaseq_tools.QualityAssessmentObject import QualityAssessmentObject
from rnaseq_tools.DatabaseObject import DatabaseObject
from rnaseq_tools import utils

def main(argv):
    # parse cmd line arguments
    print('...parsing cmd line input')
    args = parseArgs(argv)  # TODO error check all input
    try:
        if not os.path.isdir(args.reports_dir):
            raise NotADirectoryError('OutputDirDoesNotExist')
    except NotADirectoryError:
        print('%s does not lead to a valid directory. Check the path and resubmit with working -r' % args.reports)
    else:
        align_count_path = args.reports_dir
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
        run_number = utils.getRunNumber(args.reports_dir)
    except AttributeError:
        run_number = input(
            'No run number detected in input directory name. Enter something to insert in the output directory\n'
            'name: run_<what_you_input>_summary.csv: ')
    # create filename
    quality_assessment_filename = "run_{}_quality_summary.csv".format(run_number)
    output_path = os.path.join(output_directory, quality_assessment_filename)
    # create QualityAssessmentObject
    qa = QualityAssessmentObject(align_count_path=align_count_path,
                                 run_number=run_number,
                                 output_dir=args.output_dir,
                                 quality_assessment_filename=quality_assessment_filename)

    print('...compiling alignment information')
    # create dataframes storing the relevant alignment and count metadata from the novoalign and htseq logs
    qual_assess_1_df = QualityAssessmentObject.compileData(qa.align_count_path, ["_novoalign.log", "_read_count.tsv"])

    # re_order columns
    column_order = ['LIBRARY_SIZE', 'TOTAL_ALIGNMENT', 'UNIQUE_ALIGNMENT', 'MULTI_MAP', 'NO_MAP', 'HOMOPOLY_FILTER',
                    'READ_LENGTH_FILTER', 'WITH_FEATURE_RATIO', 'WITH_FEATURE', 'NO_FEATURE',
                    'NOT_ALIGNED_TO_FEATURE', 'FEATURE_ALIGN_NOT_UNIQUE', 'AMBIGUOUS_FEATURE', 'TOO_LOW_AQUAL']
    qual_assess_1_df = qual_assess_1_df[column_order]
    print('writing output to %s' % output_path)
    # write to csv
    qual_assess_1_df.to_csv(output_path, index_label="FASTQFILENAME")

    # TODO: make genotype check automatic
    # prompt user to enter gene_list if genotype check fails
    # set current to 0 for ko, <50 percentile for _over
    # make labels both genotype and fastqfilename in R script
    # take browser shot -- build this as class (certain attributes necessary to taking browser shot)

    # perturbed genotype check
    if args.coverage_check:
        print('...extracting perturbed samples from query sheet for coverage check sbatch script')
        # check if all necessary components are present
        if os.path.isfile(args.query_sheet_path):
            qa.query_sheet_path = args.query_sheet_path
            qa.query_df = utils.readInDataframe(qa.query_sheet_path)
        else:
            raise FileNotFoundError('QuerySheetDoesNotExist')
        # extract perturbed samples' COUNTFILENAME
        qa.standardized_database_df = DatabaseObject.standardizeDatabaseDataframe(qa.query_df)
        # create filter (boolean column, used in following line)
        df_wt_filter = qa.standardized_database_df['GENOTYPE'] != 'CNAG_00000'
        perturbed_sample_list = list(qa.standardized_database_df[df_wt_filter]['COUNTFILENAME'])
        # create path to new sbatch script
        sbatch_job_script_path = os.path.join(qa.job_scripts,'coverage_%s_%s.sbatch'
                                              %(qa.year_month_day, utils.hourMinuteSecond()))
        # write sbatch script
        print('...writing coverage check sbatch script')
        with open(sbatch_job_script_path, 'w') as sbatch_file:
            sbatch_file.write("#!/bin/bash\n")
            sbatch_file.write("#SBATCH --mem=5G\n")
            sbatch_file.write("#SBATCH -D %s\n" % qa.user_rnaseq_pipeline_directory)
            sbatch_file.write("#SBATCH -o sbatch_log/coverage_calculation_%A_%a.out\n")
            sbatch_file.write("#SBATCH -e sbatch_log/coverage_calculation_%A_%a.err\n")
            sbatch_file.write("#SBATCH -J coverage_calculation\n\n")
            sbatch_file.write("ml bedtools\n\n")
            for sample in perturbed_sample_list:
                sorted_alignment_file = sample.replace('_read_count.tsv', '_sorted_aligned_reads.bam')
                sorted_alignment_path = os.path.join(qa.align_count_path, sorted_alignment_file)
                coverage_filename = sample.replace('_read_count.tsv', '_coverage.tsv')
                coverage_output_path = os.path.join(qa.align_count_path, coverage_filename)
                sbatch_file.write(
                    'bedtools genomecov -ibam %s -bga > %s\n' % (sorted_alignment_path, coverage_output_path))
        print('sbatch script to quantify per base coverage in perturbed samples at %s' % sbatch_job_script_path)
        print('submitting sbatch job. Once this completes, use script quantify_perturbed_coverage.py')
        cmd = 'sbatch %s' % sbatch_job_script_path
        utils.executeSubProcess(cmd)


def parseArgs(argv):
    parser = argparse.ArgumentParser(description="This script summarizes the output from pipeline wrapper.")
    parser.add_argument("-r", "--reports_dir", required=True,
                        help="[REQUIRED] Directory for novo and htseq log files. This currently only works for Novo and htseq output.")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="[REQUIRED] Strongly suggested Usage: in reports/run_####/{organism}_pipeline_info. Remember that runs with multiple organisms will have different pipeline_info dirs per organism."
                             " File path to the directory you wish to deposit the summary. Note: the summary will be called run_###_summary.csv")
    parser.add_argument("-cc", "--coverage_check", action='store_true',
                        help="[OPTIONAL] For Crypto experiments. Set this flag to add a column to the output dataframe with percent gene coverage")
    parser.add_argument("-qs", "--query_sheet_path",
                        help="[OPTIONAL] But required with -cc is set. Path to query sheet filtered for the files contained in the path passed to -r")
    args = parser.parse_args(argv[1:])
    return args


if __name__ == "__main__":
    main(sys.argv)
