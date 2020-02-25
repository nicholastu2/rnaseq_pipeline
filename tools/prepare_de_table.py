#!/usr/bin/env python
import sys
import argparse
import os.path
import pandas as pd
import glob
import copy
import utils
from utils import *

def main(argv):
    parsed = parse_args(argv)
    if not os.path.exists(parsed.sample_summary):
        sys.exit('ERROR: %s does not exist.' % parsed.sample_summary)

    ## get conditions
    experimental_conditions = [c.strip() for c in parsed.experimental_conditions.split(',')]
    contrast_condition = [parsed.contrast_conditions]

    ## load sample summary
    summary_df = pd.read_excel(parsed.sample_summary)
    ## sift data based on group and quality
    summary_df = summary_df[(summary_df['MANUAL_AUDIT']==0)]
    ## prepare design table
    design_df = build_design_table(summary_df, experimental_conditions, contrast_condition, parsed.control_value)

    experiment_name = fileBaseName(os.path.basename(parsed.sample_summary))
    experiment_name = experiment_name.replace('_quality_summary_2', '')
    design_table_filepath = os.path.join(parsed.output_dir, experiment_name + '_design_table')
    save_dataframe(design_table_filepath, design_df)

def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sample_summary', required=True,
                        help='[Required] <your_experiment>_quality_summary_2.xlsx -- the output of quality_assess_2.py')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='[Required] Where to deposit the design table for DE anlaysis. Suggested usage: the experiment directory. The name of the file will be the basename of your sample_summary + _design_table.xlsx>')
    parser.add_argument('-e', '--experimental_conditions', default='GENOTYPE,TREATMENT',
                        help='Experimental conditions to describe the sample. Use delimiter "," if multiple descriptors are used. Default is GENOTYPE,TREATMENT')
    parser.add_argument('-c', '--contrast_conditions', default='TIMEPOINT',
                        help='Experimental conditions to describe the sample. Use delimiter "," if multiple descriptors are used. Default is currently TIMEPOINT')
    parser.add_argument('-cv',  '--control_value', default=-1,
                        help='Control, e.g. CNAG_00000 for crypto, -1 if a time course experiment. Default is currently -1')
    return parser.parse_args(argv[1:])

def build_design_table(summary_df, experimental_conditions, contrast_condition, control, timecourse=True):
    """
    Automatically build design table that has basic comparison groups.
    """
    if timecourse == True:
        experiment_design = CreateDesignMatrixColumns(experimental_conditions, contrast_condition, summary_df, control)
        return experiment_design.design_df
    else:
        pass

    # print('... Building design table')
    # ## intialize design dataframe
    # design_df = summary_df[['GENOTYPE', 'REPLICATE', 'FASTQFILENAME'] + cmp_cols].copy()
    # #cmp_cols = ['GENOTYPE'] + cmp_cols
    # ## create a dict for col: vals
    # col_dict = {}
    # for col in cmp_cols:
    # 	col_dict[col] = pd.unique(summary_df[col])
    # ## iterate thru columns to compare
    # #for col in cmp_cols:
    # for col in col_dict:
    # 	if len(col_dict[col]) < 2:
    # 		## nothing to contrast with
    # 		continue
    # 	## get value combination of other column descriptors
    # 	other_cols = copy.copy(cmp_cols)
    # 	other_cols.remove(col)
    # 	n_other_cols = len(other_cols)
    # 	other_vals = [col_dict[ocol] for ocol in other_cols]
    # 	other_vals_combos = utils.makeListProduct(other_vals)
    # 	## iterate thru with other column value fixed
    # 	for vcombo in other_vals_combos:
    # 		## get rows matching the column value combo
    # 		df2 = summary_df[sum([summary_df[other_cols[k]] == vcombo[k] \
    # 							for k in range(n_other_cols)]) == n_other_cols]
    # 		# create control statement
    # 		vcombo_name = '-'.join([':'.join([str(other_cols[k]), str(vcombo[k])]) \
    # 								for k in range(n_other_cols)])
    # 		## assess diff column types
    # 		if col == 'GENOTYPE' and control is not None:
    # 			vals = col_dict[col].tolist()
    # 			vals.remove(control)
    # 			for val in vals:
    # 				## concatenate new column and set flag
    # 				new_col = ''.join(['[',vcombo_name,']',col,':',control,'-',val]) # change delimiter btwn cntrl/value from - to &
    # 				design_df = pd.concat([design_df, pd.Series(['']*design_df.shape[0], name=new_col)], axis=1)
    # 				design_df.loc[df2.index[df2[col] == control], new_col] = '0'
    # 				design_df.loc[df2.index[df2[col] == val], new_col] = '1'
    # 		else:
    # 			vals = col_dict[col].tolist()
    # 			for val in vals:
    # 				## concatenate new column and set flag
    # 				new_col = ''.join(['[', vcombo_name, ']', col, ':', '-', str(val)])
    # 				design_df = pd.concat([design_df, pd.Series([''] * design_df.shape[0], name=new_col)], axis=1)
    # 				design_df.loc[df2.index[df2[col] == control], new_col] = '0'
    # 				design_df.loc[df2.index[df2[col] == val], new_col] = '1'
    # ## sort dataframe
    # design_df = design_df.sort_values(cmp_cols)
    # return design_df


def save_dataframe(filepath, df, df_cols=None):
    """
    Save dataframe of quality assessment.
    """
    if not filepath.endswith('.xlsx'):
        filepath += '.xlsx'
    df.to_excel(filepath, index=False, columns=df_cols)


if __name__ == '__main__':
    main(sys.argv)
