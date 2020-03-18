#!/usr/bin/env python
import sys
import os
import re
import argparse
from glob import glob
import pandas as pd
import time
import json
import utils

# current order is required -- order is hard coded into parseAlignment and parseGeneCount
ALIGN_VARS = ["READ_SEQUENCES", "UNIQUE_ALIGNMENT", "MULTI_MAPPED", "NO_MAPPING_FOUND"]
COUNT_VARS = ["TOTAL_MAPPED_READS", "WITH_FEATURE", "NO_FEATURE", "AMBIGUOUS", "TOO_LOW_A_QUAL", "NOT_ALIGNED",
              "ALIGNMENT_NOT_UNIQUE"]


def main(argv):
    """
    qual_assess_1 main method
    :param argv: cmd line input -ac path to align_counts; -rn run number; -o output directory
    """
    # parse cmd line arguments
    args = parseArgs(argv)
    # create StandardDataFormat object
    standard_data = utils.StandardDataFormat(align_count_path=args.align_count_path, run_number=args.run_number,
                                             output_dir=args.output, query_sheet_path=args.query_sheet_path,
                                             log2_cpm_path=args.log2cpm, experiment_columns=args.exp_columns)

    # If -gc not passed, this will be empty. Otherwise, this will store True
    genotype_check = args.genotype_check

    # create path to new quality_assessment_sheet
    quality_assessment_filename = "run_{}_quality_summary.csv".format(standard_data.run_number)
    output_csv = os.path.join(standard_data.output_dir, quality_assessment_filename)

    # create dataframes storing the relevant alignment and count metadata from the novoalign and htseq logs
    alignment_files_df = compileData(standard_data.align_count_path, "_novoalign.log")
    count_files_df = compileData(standard_data.align_count_path, "_read_count.tsv")

    # concat the alignment and count data
    combined_df = pd.concat([alignment_files_df, count_files_df], axis=1, sort=True, join="inner")
    # write to csv
    combined_df.to_csv(output_csv, columns=ALIGN_VARS + COUNT_VARS[1:], index_label="Sample")

    # TODO: make genotype check automatic
    # prompt user to enter gene_list if genotype check fails
    # set current to 0 for ko, <50 percentile for _over
    # make labels both genotype and fastqfilename in R script
    # take browser shot -- build this as class (certain attributes necessary to taking browser shot)

    # genotype check
    if genotype_check:
        genotypeCheck(standard_data)


def parseArgs(argv):
    parser = argparse.ArgumentParser(description="This script summarizes the output from pipeline wrapper.")
    parser.add_argument("-ac", "--align_count_path", required=True,
                        help="Directory for alignment log files. This currently only works for Novoalign output.")
    parser.add_argument("-rn", "--run_number", required=True,
                        help='the run number corresponding to this batch of fastq files')
    parser.add_argument("-o", "--output", required=True,
                        help="Suggested Usage: in reports/run_####/pipeline_info. Remember that runs with multiple organisms will have different pipeline_info dirs per organism."
                             " File path to the directory you wish to deposit the summary. Note: the summary will be called run_###_summary.csv")
    parser.add_argument("-gc", "--genotype_check",
                        help="path to a text file with a list (see templates/genotype_check.txt for an example) of fastqFileNames to perform a genotype check on")
    parser.add_argument("-qs", "--query_sheet",
                        help="path to query sheet with (only) the samples you wish to genotype check")
    parser.add_argument("-log2cpm", "--log2_cpm",
                        help="log2 cpm of the raw counts (run raw_counts.py and then log2_cpm.R to produce these")
    parser.add_argument("-exp_cols", "--experiment_columns", nargs='+',
                        help="columns to use to create replicate groups. Eg for crypto -exp_cols genotype treatment timepoint"
                             "note: no quotes or commas")
    args = parser.parse_args(argv[1:])
    return args


def compileData(dir_path, suffix):
    df = pd.DataFrame()
    file_paths = glob("{}/*{}".format(dir_path, suffix))
    for file_path in file_paths:
        fastq_filename = re.findall(r'(.+?){0}'.format(suffix), os.path.basename(file_path))[0]
        if "novoalign" in suffix:
            data = parseAlignmentLog(file_path)
        elif "read_count" in suffix:
            data = parseGeneCount(file_path)
        data.update({"Sample": fastq_filename})
        df = df.append(pd.Series(data), ignore_index=True)
    return df.set_index("Sample")


def parseAlignmentLog(alignment_log_file_path):
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
    with open(file_path, "r") as file:
        out = {COUNT_VARS[1]: 0}
        for line in file:
            k, v = line.strip().split()
            if line.startswith("__"):
                out[k[2:]] = int(v)
            else:
                out[COUNT_VARS[1]] += int(v)
        out[COUNT_VARS[0]] = sum([out[k] for k in COUNT_VARS[1:]])
        for k in COUNT_VARS[1:]:
            out[k] /= float(out[COUNT_VARS[0]])
    return out


def genotypeCheck(standard_data):
    """
    function to perform genotype check
    :param standard_data: an instance of StandardDataFormat with attributes query_sheet_path, raw_count_path, log2_cpm,
                          output_dir, experiment_columns
    :returns: None. However, outputs in subdir of align_counts a series of IGV screenshots and histograms corresponding to the genes in the
    """

    # verify standard_data has correct attributes and the paths to the query_sheet and raw_counts exit
    if not (hasattr(standard_data, 'query_sheet_path') and hasattr(standard_data, 'raw_count_path') and
            hasattr(standard_data, 'log2_cpm_path') and standard_data.experiment_columns and
            os.path.exists(standard_data.query_sheet_path) and
            os.path.exists(standard_data.raw_count_path) and os.path.exists(standard_data.log2_cpm_path)):
        sys.exit('query_sheet_path or raw_count data was either entered to StandardDataFormat incorrectly'
                 'or the paths do not exist. Check both the code constructing standard_data and the paths')
    else:
        # create script command for genotype_check_histogram.R
        script_cmd = 'genotype_check_histogram.R -q {} -c {} -o {}'.format(standard_data.standardized_query_path,
                                                                           standard_data.log2_cpm_path,
                                                                           standard_data.output_dir)
        exp_column_statement = ' -e '
        for i in standard_data.experiment_columns:
            exp_column_statement = exp_column_statement + '{},'.format(i)
        exp_column_statement = exp_column_statement[:-1]
        script_cmd = script_cmd + exp_column_statement
        # execute genotype_check_histogram.R
        os.system(script_cmd)


if __name__ == "__main__":
    main(sys.argv)
