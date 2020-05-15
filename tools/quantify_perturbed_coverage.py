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
        if not os.path.isfile(args.quality_assessment_path):
            raise FileNotFoundError('QualityAssessSheetPathNotValid')
        quality_assessment_path = args.quality_assessment_path
        quality_assessment_df = utils.readInDataframe(quality_assessment_path)
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
    od = OrganismData(organism=str(args.organism))

    df_wt_filter = standardized_query_df['GENOTYPE'] != 'CNAG_00000'
    perturbed_genotype_list = standardized_query_df[df_wt_filter]['GENOTYPE'].unique()
    # store bedfiles as dictionary
    bedfile_path_dict = {}
    for perturbed_genotype in perturbed_genotype_list:
        # put bedfiles into parent dir of quality_assess --> this is intended to be <OrganismData>_pipeline_info generally
        parent_dir_qual_assess = os.path.split(quality_assessment_path)[0]
        try:
            if not os.path.isdir(parent_dir_qual_assess):
                raise NotADirectoryError('ParentDirectoryOfQualityAssessNotFound')
        except NotADirectoryError:
            print('Could not extract the parent directory of the quality_assessment file. Try passing in the absolute path.')
        finally:
            bedfile_path = os.path.join(parent_dir_qual_assess, '%s.bed' % perturbed_genotype)
            print('checking for, and creating if DNE, %s bedfile in %s' %(perturbed_genotype, parent_dir_qual_assess))
            if not os.path.isfile(bedfile_path):
                print('...creating bed file for %s' %perturbed_genotype)
                cmd = '/home/chase/code/brentlab/rnaseq_pipeline/tools/make_bed.py -a %s -g %s -b %s' % (od.annotation_file, perturbed_genotype, bedfile_path)
                utils.executeSubProcess(cmd)
            # set key: value in dictionary
            bedfile_path_dict.setdefault(perturbed_genotype, bedfile_path)

    # create dictionary to store {perturbed_genotype: [list of coverage.tsv files]
    perturbed_coverage_dict = {}
    genotype_countfilename_df = standardized_query_df[df_wt_filter][['GENOTYPE', 'COUNTFILENAME']]  # take only these two columns
    # extract coverage files in reports_directory to match against those extracted from query_df
    perturbed_coverage_list = glob(reports_directory + '/*' + '_coverage.tsv')
    perturbed_coverage_list = [os.path.basename(coverage_file) for coverage_file in perturbed_coverage_list]

    # create dictionary to store percent coverage
    percent_coverage_dict = {}
    for index, row in genotype_countfilename_df.iterrows():
        coverage_filename = str(row['COUNTFILENAME']).replace('_read_count.tsv', '_coverage.tsv')
        if coverage_filename not in perturbed_coverage_list:
            sys.exit('coverage file extracted from query not found in the reports_directory supplied. '
                     'Please make sure that the _coverage.tsv are there, and that the query .csv is accurate to the contents of the report_directory')
        perturbed_coverage_dict.setdefault(row['GENOTYPE'], []).append(coverage_filename)

    for perturbed_genotype, coverage_path_list in perturbed_coverage_dict.items():
        for coverage_filename in coverage_path_list:
            coverage_path = os.path.join(reports_directory, coverage_filename)
            # create df from coverage file associated with sample
            per_base_coverage_df = pd.read_csv(coverage_path, sep='\t',
                                               names=['chr', 'start', 'stop', 'count'],
                                               dtype={'start': int, 'stop': int, 'count': int})
            # add column covered region (if a read spans a region in the alignment, it is counted over the region, not per base). coverage is equal to 1 * length of region covered (depth is not considered)
            per_base_coverage_df['covered_region'] = per_base_coverage_df[['count']].sum(axis=1) > 0
            per_base_coverage_df['covered_region'] = (per_base_coverage_df['stop'] - per_base_coverage_df['start']) * per_base_coverage_df['covered_region']
            # read in bed as a dataframe
            region_bed_df = pd.read_csv(bedfile_path_dict[perturbed_genotype], sep='\t', names = ['chr', 'region_start', 'region_stop', 'feature_type'])
            region_bed_df['length'] = region_bed_df.region_stop - region_bed_df.region_start
            # calculate coverage over CDS
            cds_region_df = region_bed_df[region_bed_df.feature_type == 'CDS']
            cds_region_df.reset_index(inplace=True)

            for index, row in cds_region_df.iterrows():
                cds_start = row['region_start']
                cds_stop = row['region_stop']
                chromosome = row['chr']

                cds_filter = (per_base_coverage_df['start'] >= cds_start) & (per_base_coverage_df['stop'] <= cds_stop) & (per_base_coverage_df['chr'] == chromosome)
                cds_sum = sum(per_base_coverage_df[cds_filter]['covered_region'])
                coverage = cds_sum / float(cds_region_df.loc[index, 'length'])
                cds_region_df.loc[index, 'coverage'] = coverage
            countfilename = os.path.basename(coverage_path.replace('_coverage.tsv', '_read_count.tsv'))
            percent_coverage_dict.setdefault(countfilename, sum(cds_region_df['length'] * cds_region_df['coverage']) / sum(cds_region_df['length']))

    # this is a quick possibly dirty way of determining if the quality assessment is quality_assess_1 or 2. Area of improvement.
    if 'Samples' in list(quality_assessment_df.columns):
        coverage_df = pd.DataFrame.from_dict(percent_coverage_dict, orient='index', columns=['percent_coverage']).reset_index()
        coverage_df.rename(columns={'index': 'Sample'}, inplace=True)
        coverage_df['Sample'] = coverage_df.Sample.str.replace('_read_count.tsv', '')
        qual_assess_with_percent_coverage = pd.merge(quality_assessment_df, coverage_df, on='Sample', how='left')
    else:
        coverage_df = pd.DataFrame.from_dict(percent_coverage_dict, orient='index', columns=['percent_coverage']).reset_index()
        coverage_df.rename(columns={'index': 'COUNTFILENAME'}, inplace=True)
        qual_assess_with_percent_coverage = pd.merge(quality_assessment_df, coverage_df, on='COUNTFILENAME', how='left')

    qual_assess_with_percent_coverage.to_csv(quality_assessment_path, index=False)

def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="add a column to the quality_assess_1 output quantifying gene coverage of perturbed gene")
    parser.add_argument("-r", "--reports", required=True,
                        help="[REPORTS] path to the directory with the alignment/count output if qual_assess_1 or experiment directory if qual_assess_2")
    parser.add_argument("-qa", "--quality_assessment_path", required=True,
                        help='[REQUIRED] output of quality_assess_1 or quality_assess_2.')
    parser.add_argument("-qs", "--query_sheet_path", required=True,
                        help='[REQUIRED] query sheet of only the samples in the directory passed in flag -r')
    parser.add_argument('-g', '--organism', required=True,
                        help='[REQUIRED] Must be one of organisms in genome_files and configured for use with OrganismDataObject. Currently included: KN99, S288C_R64, H99')

    args = parser.parse_args(argv[1:])
    return args


if __name__ == "__main__":
    main(sys.argv)
