#!/usr/bin/env python
import sys
import os
import re
import argparse
from glob import glob
import pandas as pd
from rnaseq_tools.QualityAssessmentObject import QualityAssessmentObject
from rnaseq_tools.DatabaseObject import DatabaseObject
from rnaseq_tools import utils

# current order is required -- order is hard coded into parseAlignment and parseGeneCount
# the formatting corresponds to the formatting in the log and read_count.tsv files
ALIGN_VARS = ["Read Sequences", "Unique Alignment", "Multi Mapped", "No Mapping Found"]
COUNT_VARS = ["total_mapped_reads", "with_feature", "no_feature", "ambiguous", "too_low_aQual", "not_aligned",
              "alignment_not_unique"]


def main(argv):
    # parse cmd line arguments
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
        print('%s does not lead to a valid directory. check the path and resubmit with correct -o' %args.output_dir)
    else:
        output_directory = args.output_dir
    try:
        if not os.path.isfile(args.query_sheet):
            raise FileNotFoundError('QuerySheetDoesNotExist')
    except FileNotFoundError:
        print('%s does not lead to a valid file. Check and resubmit correct -qs' %args.query_sheet)

    # get run number
    run_number = utils.getRunNumber(args.reports_dir)
    # create filename
    quality_assessment_filename = "run_{}_quality_summary.csv".format(run_number)
    output_path = os.path.join(output_directory, quality_assessment_filename)
    # create QualityAssessmentObject
    qa = QualityAssessmentObject(align_count_path=align_count_path,
                                 run_number=utils.getRunNumber(args.reports),
                                 output_dir=args.output,
                                 quality_assessment_filename=quality_assessment_filename, config_file='/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline_config.ini')

    # create dataframes storing the relevant alignment and count metadata from the novoalign and htseq logs
    alignment_files_df = compileData(qa.align_count_path, "_novoalign.log")
    count_files_df = compileData(qa.align_count_path, "_read_count.tsv")

    # concat the alignment and count data
    combined_df = pd.concat([alignment_files_df, count_files_df], axis=1, sort=True, join="inner")
    # write to csv
    combined_df.to_csv(output_path, columns=ALIGN_VARS + COUNT_VARS[1:], index_label="Sample")

    # TODO: make genotype check automatic
    # prompt user to enter gene_list if genotype check fails
    # set current to 0 for ko, <50 percentile for _over
    # make labels both genotype and fastqfilename in R script
    # take browser shot -- build this as class (certain attributes necessary to taking browser shot)

    # perturbed genotype check
    if args.coverage_check:
        # check if all necessary components are present
        if os.path.isfile(args.query_sheet_path):
            qa.query_sheet_path = args.query_sheet_path
        else:
            raise FileNotFoundError('QuerySheetDoesNotExist')
        # extract perturbed samples' COUNTFILENAME
        qa.standardized_database_df = DatabaseObject.standardizeDatabaseDataframe(qa.query_sheet_path)
        # create filter (boolean column, used in following line)
        df_filter = qa.standardized_database_df['GENOTYPE'] != 'CNAG_00000'
        # filter out CNAG_00000 (wt) from standardized_data_df (only perturbed genotypes remain)
        perturbed_sample_list = qa.standardized_database_df[df_filter]
        # create path to new sbatch script
        sbatch_job_script_path = os.path.join(qa.job_scripts, 'coverage_%s_%s.sbatch' %(qa.year_month_day, utils.hourMinuteSecond()))
        # write sbatch script
        with open(sbatch_job_script_path, 'w') as sbatch_file:
            sbatch_file.write("#!/bin/bash\n")
            sbatch_file.write("#SBATCH --mem=5G\n")
            sbatch_file.write("#SBATCH -D ./\n")
            sbatch_file.write("#SBATCH -o sbatch_log/coverage_calculation_%A_%a.out\n")
            sbatch_file.write("#SBATCH -e sbatch_log/coverage_calculation_%A_%a.err\n")
            sbatch_file.write("#SBATCH -J coverage_calculation\n\n")
            sbatch_file.write("ml bedtools\n")
            for sample in perturbed_sample_list:
                sorted_alignment_file = sample.replace('_read_count.tsv', '_sorted_aligned_reads.bam')
                coverage_filename = sample.replace('_read_count.tsv', '_coverage.tsv')
                sbatch_file.write('bedtools genomecov -ibam %s -bga > %s\n' %(sorted_alignment_file, coverage_filename))

def parseArgs(argv):
    parser = argparse.ArgumentParser(description="This script summarizes the output from pipeline wrapper.")
    parser.add_argument("-r", "--reports_dir", required=True,
                        help="[REQUIRED] Directory for alignment log files. This currently only works for Novoalign output.")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="[REQUIRED] Suggested Usage: in reports/run_####/{organism}_pipeline_info. Remember that runs with multiple organisms will have different pipeline_info dirs per organism."
                             " File path to the directory you wish to deposit the summary. Note: the summary will be called run_###_summary.csv")
    parser.add_argument("-cc", "--coverage_check", action='store_true',
                        help="[OPTIONAL] For Crypto experiments. Set this flag to add a column to the output dataframe with percent gene coverage")
    parser.add_argument("-qs", "--query_sheet",
                        help="[OPTIONAL] But required with -cc is set. Path to query sheet containing AT LEAST the samples you wish to genotype check. There may be more entries.")
    args = parser.parse_args(argv[1:])
    return args


def compileData(dir_path, suffix): # make this into object by itself?
    """
    get a list of the filenames in the run_#### file that correspond to a given type
    :param dir_path: path to the run_#### directory, generally (and intended to be) in /scratch/mblab/$USER/rnaseq_pipeline/reports
    :param suffix: the type of file either novoalign or _read_count
    :returns: a dataframe containing the files according to their suffix
    """
    # instantiate dataframe
    df = pd.DataFrame()
    # extract files in directory with given suffix
    file_paths = glob("{}/*{}".format(dir_path, suffix))
    for file_path in file_paths:
        # extract fastq filename
        fastq_filename = re.findall(r'(.+?)%s' %suffix, os.path.basename(file_path))[0]
        if "novoalign" in suffix:
            data = parseAlignmentLog(file_path)
        elif "read_count" in suffix:
            data = parseGeneCount(file_path)
        data.update({"Sample": fastq_filename})
        df = df.append(pd.Series(data), ignore_index=True)
    return df.set_index("Sample")


def parseAlignmentLog(alignment_log_file_path):
    """
    parse the valuable informatino out of the alignment logs
    :param alignment_log_file_path: the filepath to a novoalign alignment log
    :returns: a dictionary of the parsed data of the input file
    """
    with open(alignment_log_file_path, "r") as file:
        alignment_metadata_dict = {}
        for line in file:
            line = line.strip("#").strip()
            for k in ALIGN_VARS:
                if line.startswith(k):
                    v = int(line.split(":")[-1].split("(")[0].strip())
                    alignment_metadata_dict[k] = v if k == ALIGN_VARS[0] else v / float(
                        alignment_metadata_dict[ALIGN_VARS[0]])
    return alignment_metadata_dict


def parseGeneCount(file_path):
    """
    count the gene counts that mapped either to genes or to other features (see COUNT_VARS at top of script for other features)
    :param file_path: a path toa  _read_count.tsv file (htseq output)
    :returns: a dictionary with the features in COUNT_VARS quantified for the inputted file
    """
    with open(file_path, "r") as file:
        # COUNT_VARS[1] is WITH_FEATURE
        out = {COUNT_VARS[1]: 0}
        # walk through the count file
        for line in file:
            # split the two columns and store values as index and value
            index, value = line.strip().split()
            # if the line starts with __, it stores information such as __no_feature
            if line.startswith("__"):
                out[index[2:]] = int(value)
            else:
                # otherwise the line will store individual gene counts. Add these counts to WITH_FEATURE
                out[COUNT_VARS[1]] += int(value)
        # COUNT_VARS[0] = TOTAL_MAPPED_READS and is the sum across the rows of the counts in the various categories
        out[COUNT_VARS[0]] = sum([out[key] for key in COUNT_VARS[1:]])
        # divide each category starting with WITH_FEATURE by the TOTAL_MAPPED_READS
        for key in COUNT_VARS[1:]:
            out[key] /= float(out[COUNT_VARS[0]])
    return out


if __name__ == "__main__":
    main(sys.argv)
