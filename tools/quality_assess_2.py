#!/usr/bin/env python
import sys
import argparse
import re
import pandas as pd
import numpy as np
from rnaseq_tools.OrganismDataObject import OrganismData
from rnaseq_tools import utils
import os


def main(argv):
    # parse cmd line input and instantiate OrganismData
    parsed = parse_args(argv)
    od = OrganismData(organism=parsed.organism, query_sheet_path=parsed.query_sheet_path,
                      experiment_dir=parsed.experiment_dir, norm_count_path=parsed.norm_count_path,
                      max_replicates=parsed.max_replicates, output_dir=parsed.output_dir,
                      wildtype=parsed.wildtype, drug_marker=parsed.drug_marker, qc_config=parsed.qc_config,
                      experiment_conditions=parsed.experimental_conditions.split(' '), interactive=True)  # TODO -- deal with multiple inputs better than this. point of weakness
    # TODO: put in flag to check if not on HTCF
    od.standardDirectoryStructure()
    od.createOrganismDataLogger()
    # create output sheet name
    experiment_dir = utils.dirName(od.experiment_dir)
    filename = experiment_dir + '_quality_summary_2.xlsx'
    output_name = os.path.join(od.output_dir, filename)

    # validate paths
    if os.path.exists(output_name):
        print('WARNING: %s already exists. If you continue, the file will be overwritten. '
              'Do you wish to continue (y/n)?\n' % output_name)
        user_response = input()
        if user_response == 'n':
            sys.exit('goodbye')
    if not os.path.exists(od.norm_count_path):
        sys.exit('ERROR: %s (the normalized count matrix) does not exist.' % od.norm_count_path)
    if not os.path.exists(os.path.dirname(output_name)):
        sys.exit('ERROR: %s (the output directory) does not exist.' % os.path.dirname(output_name))

    # load QC config data
    # TODO: complexity.thresh <- mean(alignment.sum$COMPLEXITY[indx]) - 2*sd(alignment.sum$COMPLEXITY[indx]);
    global QC_dict  # TODO: change this to load directly to od object
    QC_dict = utils.loadConfig(od.qc_config)

    # # get conditions #
    # conditions = None if parsed.condition_descriptors is None else \  # this should be taken care of by nargs = '+' -- check
    #     [c.strip() for c in parsed.condition_descriptors.split(',')]

    ## do QA
    print('... Preparing QA dataframe')
    if od.drug_marker is None:
        drug_markers = None
        drug_marker_columns = []
    else:
        drug_markers = od.drug_marker
        drug_marker_columns = [drug_marker + '_FOM' for drug_marker in od.drug_marker]
    df_columns = ['GENOTYPE', 'REPLICATE', 'COUNTFILENAME'] \
                 + od.experiment_conditions \
                 + ['STATUS', 'AUTO_AUDIT', 'MANUAL_AUDIT', 'USER', 'NOTE'] \
                 + ['TOTAL', 'ALIGN_PCT', 'MUT_FOW'] \
                 + drug_marker_columns \
                 + ['COV_MED_REP' + ''.join(np.array(combo, dtype=str)) for combo
                    in utils.makeCombinations(range(1, od.max_replicates + 1))]
    qual_assess_df, rep_max = initializeQualAssesDf(od.query_df, df_columns, od.experiment_conditions)
    if rep_max != od.max_replicates:
        print(
            'The max number of replicates in the query sheet is {}. Please re-launch this script with -r {}. However, calculating CoV with greater than 7 samples is not possible currently.'.format(
                rep_max, rep_max))
    if rep_max > 7:
        print(
            'Calculating the power set of 7 replicates for the CoV is not currently possible as it exceeds the pandas maximum sheet size. A QA sample summary will be returned, but CoV will not be calculated. Do you wish to continue? Enter y or n')
        user_response = input()
        if user_response == 'n':
            sys.exit()
        else:  # TODO: clean this up -- this is repeat code to avoid assess_replicate_concorndance if number of samples too large
            expr, sample_dict = loadExpressionData(qual_assess_df, od.norm_count_path, od.gene_list,
                                                    od.experiment_conditions)
            print('... Assessing reads mapping')
            df = assess_mapping_quality(qual_assess_df, od.experiment_dir)
            print('... Assessing efficiency of gene mutation')
            if parsed.descriptors_specific_fow:
                df = assess_efficient_mutation(df, expr, sample_dict, od.wildtype,
                                               od.experiment_conditions)  # TODO: clean up the df's in here -- probably should be qual_assess_df
            else:
                df = assess_efficient_mutation(df, expr, sample_dict, od.wildtype)
            print('... Assessing insertion of resistance cassette')
            df = assess_resistance_cassettes(df, expr, drug_markers, od.wildtype)
            df = update_auto_audit(df, parsed.auto_audit_threshold)
            save_dataframe(output_name, df, df_columns, od.experiment_conditions, len(od.experiment_conditions))
            sys.exit()
    norm_count_df, sample_dict = loadExpressionData(qual_assess_df, od.norm_count_path, od.gene_list,
                                                    od.experiment_conditions)
    print('... Assessing reads mapping')
    qual_assess_df = assess_mapping_quality(qual_assess_df, od.output_dir)
    print('... Assessing efficiency of gene mutation')
    if parsed.descriptors_specific_fow:
        qual_assess_df = assess_efficient_mutation(qual_assess_df, norm_count_df, sample_dict, od.wildtype,
                                                   od.experiment_conditions)
    else:
        qual_assess_df = assess_efficient_mutation(qual_assess_df, norm_count_df, sample_dict, od.wildtype)
    print('... Assessing insertion of resistance cassette')
    qual_assess_df = assess_resistance_cassettes(qual_assess_df, norm_count_df, drug_markers, od.wildtype)
    print('... Assessing concordance among replicates')
    #qual_assess_df = assess_replicate_concordance(qual_assess_df, norm_count_df, sample_dict, od.experiment_conditions)
    print('... Auto auditing')
    qual_assess_df = update_auto_audit(qual_assess_df, parsed.auto_audit_threshold)
    print('...writing summary to %s' % output_name)
    save_dataframe(output_name, qual_assess_df, df_columns, od.experiment_conditions, len(od.experiment_conditions))


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-qs', '--query_sheet_path', required=True,
                        help='[REQUIRED] The output of queryDB that was used to create this experiment')
    parser.add_argument('-g', '--organism', required=True,
                        help='[REQUIRED] The organism (corresponds to the genomes_files organisms')
    parser.add_argument('-e', '--experiment_dir', required=True,
                        help='[REQUIRED] the path to the experiment directory created by create_experiment')
    parser.add_argument('-c', '--norm_count_path', required=True,
                        help='[REQUIRED] Normalized count matrix. If not given, the filepath will be guessed based on '
                             'analysis group number.')
    parser.add_argument('-r', '--max_replicates', required=True, type=int,
                        help='[REQUIRED] Maximal number of replicate in experiment design.')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='[REQUIRED] directory in which to deposit the sample quality summary.')
    parser.add_argument('-w', '--wildtype',
                        help='Wildtype genotype, e.g. CNAG_00000 for crypto, BY4741 for yeast.')
    parser.add_argument('-d', '--drug_marker', nargs='+', default=None,
                        help='drug markers inserted to replace the deleted genes. List without delimiter if multiple exist '
                             'eg CNAG_G418 CNAG_NAT')
    parser.add_argument('--experimental_conditions', nargs='+', default='TREATMENT TIMEPOINT',
                        help='[Default TREATMENT TIMEPOINT] Experimental conditions that describe the sample are used to identify subgroups within each '
                             'genotype. No delimiter necessary if multiple conditions used.')
    parser.add_argument('--descriptors_specific_fow', action='store_true',
                        help='Set this flag to find the wildtype samples that match the condition descriptors of '
                             'the mutant sample when calculating the fold change over wildtype (FOW).')
    parser.add_argument('--qc_config',
                        default='/opt/apps/labs/mblab/software/rnaseq_pipeline/1.0/config/qc_config.yaml',
                        help='Configuration file for quality assessment.')
    parser.add_argument('--auto_audit_threshold', type=int, default=0,
                        help='Threshold for automatic sample audit.')
    return parser.parse_args(argv[1:])


def initializeQualAssesDf(standardized_query_df, df_cols, conditions):
    """
    Define the qual_assess_2 dataframe
    :param standardized_query_df: A standardized query, manupulated by a StandardData object (ie OrganismData in this case).
    See StandardData for info on what standard_query format is
    :param df_cols: columns for the quality_assessment_dataframe
    :param conditions: additional conditions to append to the quality assessment from the standardized_query_sheet
    :returns: the blank (except for the query sheet columns) sample_summary and the maximum number of replicates
    """
    quality_assess_2_df = pd.DataFrame(columns=df_cols)

    # cast replicate to int
    standardized_query_df = standardized_query_df.astype({'REPLICATE': 'float'})
    standardized_query_df = standardized_query_df.astype({'REPLICATE': 'int32'})

    # if passed, force conditions to be appended to sample_summary_df from standardized_query_df to uppercase
    if conditions:
        conditions = [x.upper() for x in conditions]

    standardized_query_df = standardized_query_df[['GENOTYPE', 'REPLICATE', 'COUNTFILENAME'] + conditions]
    standardized_query_df = standardized_query_df.reset_index().drop(['index'], axis=1)
    standardized_query_df = pd.concat(
        [standardized_query_df, pd.Series([0] * standardized_query_df.shape[0], name='STATUS')], axis=1)
    standardized_query_df = pd.concat(
        [standardized_query_df, pd.Series([np.nan] * standardized_query_df.shape[0], name='AUTO_AUDIT')], axis=1)

    # create sample_summary dataframe
    quality_assess_2_df = quality_assess_2_df.append(standardized_query_df)
    # re-index replicates. This is absolutely necessary b/c there may be technical replicates with the same biological replicate number
    conditions.append('GENOTYPE')
    quality_assess_2_df['REPLICATE'] = quality_assess_2_df.groupby(conditions).cumcount() + 1
    # get maximum number of replicates -- use this as a check against what the user entered
    rep_max = quality_assess_2_df.groupby(conditions).cumcount().max() + 1

    return quality_assess_2_df, rep_max


def loadExpressionData(qual_assess_df, norm_count_matrix, gene_list, experiment_conditions):
    """
    Load count matrix, and make a sample dictionary.
    :param qual_assess_df:
    :param norm_count_matrix: path to the norm_count.csv
    :param gene_list: list of gene names -- these are found in genome_files in the subdirectory of the given organism. This is a attribute of any OrganismData object
    :param experiment_conditions: control conditions of experiment
    :returns:
    """
    # load count matrix
    norm_count_df = pd.read_csv(norm_count_matrix)
    norm_count_df = norm_count_df.rename(columns={'Unnamed: 0': 'gene'})

    # intersect gene list with the count matrix, remove genes from count matrix that are not in gene list
    if gene_list is not None:
        gene_list = pd.read_csv(gene_list, names=['gene'])
        if len(np.setdiff1d(gene_list, norm_count_df['gene'])) != 0:
            print('WARNING: The custom gene list contains genes that are not in count matrix. '
                  'Proceeding using the intersection.')
        # convert to array
        gene_list = np.intersect1d(gene_list, norm_count_df['gene'])
        # remove genes NOT IN the gene list from the count matrix
        norm_count_df = norm_count_df.loc[norm_count_df['gene'].isin(gene_list)]

    # make sample dict with structure {(genotype, condition1, condition2, ...): {rep#: countfilename}} eg { ('CNAG_00883', '37C.CO2', 90): {1: 'Brent_3_GTAC_2_SIC_Index2_07_GCTTAGAA_GAGTTGGT_S4_R1_001_read_count.tsv', 2: 'Brent_4_GTAC_4_SIC_Index2_07_CACCTCCA_GAGTTGGT_S5_R1_001_read_count.tsv', 3: 'Brent_5_GTAC_5_SIC_Index2_07_ATCGAGCA_GAGTTGGT_S6_R1_001_read_count.tsv'}}
    sample_dict = {}
    for index, row in qual_assess_df.iterrows():
        genotype = row['GENOTYPE']
        # set sample description, eg (genotype, condition1, condition2, ...)
        if len(experiment_conditions) == 0:
            sample_description = tuple([genotype])
        else:
            sample_description = tuple([genotype] + [row[c] for c in experiment_conditions])
        count_file_name = str(row['COUNTFILENAME'])

        # that the count_file_name is in fact in the norm counts
        if count_file_name in norm_count_df.columns.values:
            # and then add to the dictionary
            if sample_description not in sample_dict.keys():
                sample_dict[sample_description] = {}
            sample_dict[sample_description][row['REPLICATE']] = count_file_name

    return norm_count_df, sample_dict


def assess_mapping_quality(qual_assess_df, experiment_directory, aligner_tool='novoalign'):
    """
    Assess percentage of uniquely mapped reads over all reads
    :param qual_assess_df: the work-in-progress quality_assessment_df
    :param experiment_directory: directory where the count and log files live (should be in OrganismData.experiment_dir
    :param aligner_tool: the aligner tool used. currently, only novoalign. this is the suffix.log that will be appended
    to search for log files in the experiment_dir
    :returns: updated quality_assess_df
    """
    for i, row in qual_assess_df.iterrows():
        sample = str(row['COUNTFILENAME']).replace('_read_count.tsv', '_%s.log' % aligner_tool)
        filepath = os.path.join(experiment_directory, sample)

        # read alignment log
        with open(filepath, 'r') as reader:
            lines = reader.readlines()
            for line in lines:
                reg_total = re.search(r'Read Sequences:( +)(.\d+)', line)
                reg_uniq = re.search(r'Unique Alignment:( +)(.\d+)', line)
                if reg_total:
                    total_reads = int(reg_total.group(2))
                if reg_uniq:
                    uniq_mapped_reads = int(reg_uniq.group(2))

        align_pct = uniq_mapped_reads / float(total_reads)
        # set mapping quality
        row['TOTAL'] = total_reads  # read Sequences
        row['ALIGN_PCT'] = align_pct  # Unique Alignment
        if total_reads < QC_dict['TOTAL_READS']['threshold']:
            row['STATUS'] += QC_dict['TOTAL_READS']['status']
        if align_pct < QC_dict['ALIGN_PCT']['threshold']:
            row['STATUS'] += QC_dict['ALIGN_PCT']['status']
        qual_assess_df.iloc[i] = row

    return qual_assess_df


def assess_efficient_mutation(qual_assess_df, norm_count_df, sample_dict, wt, conditions=None):
    """
    Assess the completeness of gene deletion or efficiency of gene overexpression by caluclating the expression
    of the perturbed genein mutant sample over mean expression of the same gene in wildtype.
    :param qual_assess_df: the quality_assessment_df in progress
    :param norm_count_df: normalized count dataframe
    :param sample_dict: see loadExpressionData() above for extensive description
    :param wt: wildtype ie if crypto CNAG_00000
    :param conditions: experimental conditions (ie timepoint treatment)
    :returns: updated quality_assessment_df
    """
    if wt is None:
        return qual_assess_df
    descr_match = True if conditions is not None else False
    # get wildtype samples if not matching descriptors
    if not descr_match:
        wt_samples = []
        for key in sample_dict.keys():
            if key[0] == wt:
                wt_samples += sample_dict[key].values()
        # calculate mean expression level of each gene
        wt_expr = pd.Series(pd.DataFrame.mean(norm_count_df[wt_samples], axis=1), name='mean_fpkm')
        wt_expr = pd.concat([norm_count_df, wt_expr], axis=1)
    # calculate efficiency of gene deletion, ignoring overexpression(*_over)
    for i, row in qual_assess_df[qual_assess_df['GENOTYPE'] != wt].iterrows():
        sample = str(row['COUNTFILENAME'])
        print('... calculating MUT_FOW in %s' % sample)
        ## check for each mutant gene (there could be multiple mutant genes, delimited by '.')
        mut_fow_list = []
        for mut_gene in row['GENOTYPE'].split('.'):
            # get wildtype samples if not matching descriptors
            if descr_match:
                mut_descr = [row[c] for c in conditions]
                wt_samples = []
                for key in sample_dict.keys():
                    descr_matched = all([key[j + 1] == mut_descr[j] for j in range(len(mut_descr))])
                    if key[0] == wt and descr_matched:
                        wt_samples += sample_dict[key].values()
                if len(wt_samples) == 0:
                    print(
                        '\tSample %s has no WT sample that matches its condition descriptors. Skipping this sample' % sample)
                    continue
                ## calculate mean expression level of each gene
                wt_expr = pd.Series(pd.DataFrame.mean(norm_count_df[wt_samples], axis=1), name='mean_fpkm')
                wt_expr = pd.concat([norm_count_df, wt_expr], axis=1)
            ## get mutant gene expression in mutatnt sample
            mut_gene2 = mut_gene.strip("_over")
            if mut_gene2 not in norm_count_df['gene'].tolist():
                print('\t%s not in gene list. Skipping this genotype' % mut_gene2)
                continue
            wt_mean = float(wt_expr[wt_expr['gene'] == mut_gene2]['mean_fpkm'])
            if wt_mean == 0:
                print('\t%s has 0 mean expression in WT samples' % mut_gene2)
                mut_fow = np.inf
            else:
                mut_fow = float(norm_count_df[norm_count_df['gene'] == mut_gene2][sample]) / wt_mean

            if mut_gene.endswith('_over'):
                ## check overexpression
                if (mut_fow < QC_dict['MUT_FOW']['OVEREXPRESSION']['threshold']) and (
                        row['STATUS'] < QC_dict['MUT_FOW']['OVEREXPRESSION']['status']):
                    row['STATUS'] += QC_dict['MUT_FOW']['OVEREXPRESSION']['status']
            else:
                ## check deletion
                if (mut_fow > QC_dict['MUT_FOW']['DELETION']['threshold']) and (
                        row['STATUS'] < QC_dict['MUT_FOW']['DELETION']['status']):
                    row['STATUS'] += QC_dict['MUT_FOW']['DELETION']['status']
            mut_fow_list.append(str(mut_fow))
        row['MUT_FOW'] = ','.join(mut_fow_list)
        qual_assess_df.iloc[i] = row
        print('complete. Moving onto next sample.')
    return qual_assess_df


def assess_resistance_cassettes(df, norm_count_df, resi_cass, wt):
    """
    Assess drug resistance marker gene expression, making sure the proper
    marker gene is swapped in place of the perturbed gene.
    """
    if resi_cass is None:
        return df
    ## get the median of resistance cassettes
    rc_med_dict = {}
    mut_samples = [s for s in norm_count_df.columns.values if (not s.startswith(wt)) and (s != 'gene')]
    for rc in resi_cass:
        ## exclude wildtypes and markers expressed < 150 normalized counts
        rc_fpkm = norm_count_df.loc[norm_count_df['gene'] == rc, mut_samples] # recall that the index column for norm_count_df was renamed in loadExpressionData to 'gene'
        rc_fpkm = rc_fpkm.loc[:, (np.sum(rc_fpkm, axis=0) > 150)]
        rc_med_dict[rc] = rc_fom = np.nan if rc_fpkm.empty else np.median(rc_fpkm)
    ## calcualte FOM (fold change over mutant) of the resistance cassette
    for i, row in df.iterrows():
        genotype = row['GENOTYPE']
        sample = str(row['COUNTFILENAME'])
        print('...assessing_drug_marker in %s' % genotype)
        ## update FOM
        for rc in rc_med_dict.keys():
            row[rc + '_FOM'] = np.nan if np.isnan(rc_med_dict[rc]) else float(norm_count_df.loc[norm_count_df['gene'] == rc, sample]) / \
                                                                        rc_med_dict[rc]
        ## flag those two problems:
        ## 1. the resistance cassette is exprssed in WT
        ## 2. more than one resistance cassette is expressed in single mutant
        ## TODO: add criteria for multi-mutants
        fom_check = [row[rc + '_FOM'] > QC_dict['MARKER_FOM']['threshold'] * rc_med_dict[rc] for rc in
                     rc_med_dict.keys()]
        if genotype == wt and sum(fom_check) > 0:
            row['STATUS'] += QC_dict['MARKER_FOM']['status']
        if genotype != wt and len(genotype.split('.')) > 1 and sum(fom_check) > 1:
            row['STATUS'] += QC_dict['MARKER_FOM']['status']
        print('complete. moving onto next sample.')
    return df


def assess_replicate_concordance(qual_assess_df, expr, sample_dict, conditions):
    """
    Assess the concordance among the replicates of each genotype by calculating
    the COV of each combination of replicates. Then find the maximal number of
    concordant replicates.
    """

    # calculate COV medians for replicate combinations
    for key in sorted(sample_dict.keys()):
        sample_ids = [s for s in sample_dict[key].items()]
        if len(sample_ids) == 1:
            continue # only one replicate, cannot calculate cov
        cov_meds_dict = {}
        rep_combos = utils.makeCombinations(sample_dict[key].keys())

        for rep_combo in rep_combos:
            # sort as integers
            # rep_combo = np.array(sorted(np.array(rep_combo,dtype=int)), dtype=str)
            rep_num = len(rep_combo)
            sample_combo = [sample_dict[key][rep] for rep in rep_combo]
            # calculate COV median
            cov_median = calculate_cov_median(expr[sample_combo])
            rep_combo_col = 'COV_MED_REP' + ''.join(np.array(rep_combo, dtype=str))
            qual_assess_df.loc[qual_assess_df['COUNTFILENAME'].isin(sample_combo), rep_combo_col] = cov_median
            # store COV median at the respective rep number
            if rep_num not in cov_meds_dict.keys():
                cov_meds_dict[rep_num] = {'rep_combos': [], 'cov_meds': []}
            cov_meds_dict[rep_num]['rep_combos'].append(rep_combo)
            cov_meds_dict[rep_num]['cov_meds'].append(cov_median)

        # if there are only two replicates, continue to next rep_combo in rep_combos
        if len(rep_combos) <= 1:
            continue

        ## find the maximal number of replicates that pass concordance threshold
        for rep_num in sorted(cov_meds_dict.keys())[::-1]:  # what does [:: mean?
            rep_combo = cov_meds_dict[rep_num]['rep_combos']
            cov_meds = cov_meds_dict[rep_num]['cov_meds']
            if sum([c < QC_dict['COV_MED']['threshold'] for c in cov_meds]) > 0:  # TODO: ERROR HERE
                best_combo = rep_combo[np.argmin(cov_meds)]
                break

        max_rep_combo = cov_meds_dict[max(cov_meds_dict.keys())]['rep_combos'][0]
        outlier_reps = set(max_rep_combo) - set(best_combo)

        # update STATUS column -- see qc_config. STATUS is encoded in a bit code
        for rep in outlier_reps:
            outlier_indx = set(qual_assess_df.index[(qual_assess_df['GENOTYPE'] == key[0]) & \
                                                    (qual_assess_df['REPLICATE'] == rep)])
            for ci in range(len(conditions)):
                outlier_indx = outlier_indx & \
                               set(qual_assess_df.index[qual_assess_df[conditions[ci]] == key[ci + 1]])
            qual_assess_df.loc[list(outlier_indx), 'STATUS'] += QC_dict['COV_MED']['status']

    return qual_assess_df


def calculate_cov_median(x):
    """
    Calculate the median of COVs (coefficient of variation) among replicates
    """
    covs = np.std(x, axis=1) / np.mean(x, axis=1)
    return np.nanmedian(covs)


def update_auto_audit(df, threshold):
    """
    Automatically flag sample with status over threshold
    """
    df.loc[df['STATUS'] > threshold, 'AUTO_AUDIT'] = 1
    return df


def save_dataframe(filepath, df, df_cols, conditions, fp_ext=0):
    """
    Save dataframe of quality assessment
    """
    df = df.sort_values(['GENOTYPE'] + conditions + ['REPLICATE'])
    # if not filepath.endswith('.xlsx'):
    #	filepath += '.xlsx'

    # TODO: work this in earlier in the process
    for index, row in df.iterrows():
        # df['COUNTFILENAME'] = df['COUNTFILENAME'].apply(lambda row: utils.fileBaseName(row) + '_read_count.tsv')
        pass
    df.to_excel(filepath, columns=df_cols, index=False, freeze_panes=(1, 3 + fp_ext))


if __name__ == '__main__':
    main(sys.argv)
