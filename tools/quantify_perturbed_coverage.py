#!/usr/bin/env python3
"""
   using gtf annotation file and bedtools genomecov -bga output (see coverage check in qual_asses_1), quantify
   positions in perturbed gene CDS with at least 1 read. percent_coverage = number of bases in CDS with at least 1 read / total number of bases in CDS
   usage: quantify_perturbed_coverage -qa /path/to/quality_assess_1_output.csv -g <Organism in genome files and configured for OrganismDataObject eg KN99>
"""
import sys
import argparse
import os
from glob import glob
import pandas as pd
from rnaseq_tools import utils
from rnaseq_tools.DatabaseObject import DatabaseObject
from rnaseq_tools.OrganismDataObject import OrganismData


def main(argv):
    """ main method
    :param argv: cmd line arguments
    """
    # parse cmd line arguments
    args = parseArgs(argv)
    # TODO: add error checking for organism in OrganismData object itself
    try:
        if not os.path.isfile(args.quality_assess_1_path):
            raise FileNotFoundError('QualityAssessSheetPathNotValid')
        quality_assessment_path = args.quality_assess_1_path
    except FileNotFoundError:
        sys.exit('path to quality_assess_1 not valid')
    try:
        if not os.path.isfile(args.query_sheet_path):
            raise FileNotFoundError('QuerySheetPathNotValid')
        query_sheet_path = args.query_sheet_path
        query_df = utils.readInDataframe(query_sheet_path)
        standardized_query_df = DatabaseObject.standardizeDatabaseDataframe(query_df)
    except FileNotFoundError:
        sys.exit('path to query sheet not valid')
    try:
        if not os.path.isdir(args.reports):
            raise NotADirectoryError('ReportsDirPathNotValid')
        reports_directory = args.reports
    except NotADirectoryError:
        sys.exit('path to reports directory is not valid')

    # create OrganismData object. Since organism is passed. setOrganismData is called which includes StandardData.standardDirectoryStructure(). Pass config_file to this constructor if not running on htcf
    od = OrganismData(organism=args.organism,
                      config_file='/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline_config.ini')

    print('...creating bedfiles for perturbed genotype')
    df_wt_filter = standardized_query_df['GENOTYPE'] != 'CNAG_00000'
    perturbed_genotype_list = standardized_query_df[df_wt_filter]['GENOTYPE'].unique()
    # store bedfiles as dictionary
    bedfile_path_dict = {}
    for perturbed_genotype in perturbed_genotype_list:
        # put bedfiles into parent dir of quality_assess --> this is intended to be <OrganismData>_pipeline_info generally
        parent_dir_qual_assess = utils.dirPath(quality_assessment_path)
        bedfile_path = os.path.join(parent_dir_qual_assess, '%s.bed' % perturbed_genotype)
        if not os.path.isfile(bedfile_path):
            parsed_gtf = '%s_parsed.gtf' % perturbed_genotype
            cmd = 'cat %s | grep %s > %s' % (od.annotation_file, perturbed_genotype, parsed_gtf)
            utils.executeSubProcess(cmd)


        # set key: value in dictionary
        bedfile_path_dict.setdefault(perturbed_genotype, bedfile_path)

    # create dictionary to store {perturbed_genotype: [list of coverage.tsv files]
    perturbed_coverage_dict = {}
    genotype_countfilename_df = standardized_query_df[df_wt_filter][
        ['GENOTYPE', 'COUNTFILENAME']]  # take only these two columns
    perturbed_coverage_list = glob(
        reports_directory + '/*' + '_coverage.tsv')  # extract coverage files in reports_directory to match against those extracted from query_df
    for index, row in genotype_countfilename_df.iterrows():
        coverage_filename = str(row['COUNTFILENAME']).replace('_read_count.tsv', '_coverage.tsv')
        if coverage_filename not in perturbed_coverage_list:
            sys.exit('coverage file extracted from query not found in the reports_directory supplied. '
                     'Please make sure that the _coverage.tsv are there, and that the query .csv is accurate to the contents of the report_directory')
        perturbed_coverage_dict.setdefault(row['GENOTYPE'], []).extend(coverage_filename)

    for coverage_path in perturbed_coverage_list:
        per_base_coverage_df = pd.read_csv(coverage_path, sep='\t',
                                           names=['chr', 'start', 'stop', 'count'],
                                           dtype={'start': int, 'stop': int, 'count': int})


def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="add a column to the quality_assess_1 output quantifying gene coverage of perturbed gene")
    parser.add_argument("-r", "--reports", required=True,
                        help="[REPORTS] path to the directory with the alignment and count logs where the coverage files were deposited by qual_assess_1. most likely reports/run_####")
    parser.add_argument("-qa", "--quality_assess_1_path", required=True,
                        help='[REQUIRED] output of quality_assess_1. Should be located in reports/run_####/<organism>_pipeline_info')
    parser.add_argument("-qs", "--query_sheet_path", required=True,
                        help='[REQUIRED] query sheet -- should be the same one submitted to quality_assess_1 for coverage check')
    parser.add_argument('-g', '--organism', required=True,
                        help='[REQUIRED] Must be one of organisms in genome_files and configured for use with OrganismDataObject. Currently included: KN99, S288C_R64, H99')

    args = parser.parse_args(argv[1:])
    return args


if __name__ == "__main__":
    main(sys.argv)
