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
	if not os.path.exists(parsed.samples):
		sys.exit('ERROR: %s does not exist.' % parsed.samples)

	## get conditions
	conditions = [c.strip() for c in parsed.condition_descriptors.split(',')]

	## load sample summary
	summary_df = pd.read_excel(parsed.samples)
	## sift data based on group and quality
	summary_df = summary_df[(summary_df['MANUAL_AUDIT']==0)]
	## prepare design table
	design_df = build_design_table(summary_df, conditions, parsed.wildtype)

	design_table_filepath = os.path.join(parsed.output_dir, parsed.experiment_name + '_design_table')
	save_dataframe(design_table_filepath, design_df)

def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--samples', required=True,
						help='[Required] Sample summary metadata file.')
	parser.add_argument('-g', '--group_num', required=True,
						help='[Required] Analysis group number. It is required if mutilple metadata files are used.')
	parser.add_argument('-o', '--output_dir', required=True,
						help='[Required] Auto-gen design table for DE anlaysis.')
	parser.add_argument('-n', '--experiment_name', required=True,
						help="[Required] the name of your experiment (this will be used for naming the DE table output)")
	parser.add_argument('-c',  '--control',
						help='[Required] Control, e.g. CNAG_00000 for crypto, BY4741 for yeast, EtoH if control is a condition, -1 if a timePoint.')
	parser.add_argument('--condition_descriptors', default='TREATMENT,TIMEPOINT',
						help='Experimental conditions to describe the sample. Use delimiter "," if multiple descriptors are used. Default is TREATMENT,TIME_POINT')
	return parser.parse_args(argv[1:])

# work in progress -- clean up build_design_matrix
def buildDesignTable(summary_df, contrast_cols, control_column_name, control_value):
	"""
	Build design table with basic comparison groups.
	Args: the output of quality_assessment_2.py (typically <experiment_name>_quality_summary_2.csv.
	      Note that this must be manually audited and 0s added to the manual audit column.
	      contrast_cols are the columns to be contrasted in the DE analysis
	      wt is None by default. If passed at cmd line, this will enter 0s in the design matrix across all wt rows.
	Returns: a experiment design table identifying contrast groups. See templates in the rnaseq_pipeline for an example.
	"""
	print('... Building design table')
	design_df = summary_df[['GENOTYPE', 'REPLICATE'] + contrast_cols].copy()

	# get list of columns in the contrast group
	#contrast_cols = ['GENOTYPE'] + contrast_cols

	# create a dict of the possible contract groups (eg {'GENOTYPE': ['GCN4, 'ABC1'], 'TREATMENT': ['EtoH', 'ESTRADIOL']}
	contrast_column_dict = {}
	for col in contrast_cols:
		contrast_column_dict[col] = pd.unique(summary_df[col])

	# There is a lot of repetition in this loop -- needs to be cleaned up/streamlined. point of weakness.
	for col in contrast_column_dict:
		if len(contrast_column_dict[col]) < 2:
			# nothing to contrast with
			continue # iterate to next col in contrast_cols
		else:
			cols_to_contrast_dict = contrast_column_dict.copy()
			cols_to_contrast_dict.pop(col)
			cols_to_contrast = list(cols_to_contrast_dict.keys())
			# used to be other_vals_combos
			contrast_groups = utils.makeListProduct(list(cols_to_contrast_dict.values()))
			contrast_group_count = len(cols_to_contrast_dict)
			for contrast_values in contrast_groups:
				# inside the list comprehension, return true in rows where the contrast_value matches the dataframe.
				# If the sum of the booleans in the row equal the countrast_group_count
				temp_design_matrix = summary_df[sum([summary_df[cols_to_contrast[k]] == contrast_values[k] \
													 for k in range(contrast_group_count)]) == contrast_group_count]
				# create new column headers for all contrast groups
				vcombo_name = '-'.join([':'.join([str(cols_to_contrast[k]), str(contrast_values[k])]) for k in range(contrast_group_count)])
			## assess diff column types
			if col == control_column_name and control_value is not None:
				vals = cols_to_contrast[col].tolist()
				vals.remove(control_value)
				for val in vals:
					## concatenate new column and set flag
					new_col = new_col = ''.join(['[',vcombo_name,']',col,':',control_value,'-',val])
					design_df = pd.concat([design_df, pd.Series(['']*design_df.shape[0], name=new_col)], axis=1)
					design_df.loc[temp_design_matrix.index[temp_design_matrix[col] == control_value], new_col] = '0'
					design_df.loc[temp_design_matrix.index[temp_design_matrix[col] == val], new_col] = '1'
			else:
				vals = cols_to_contrast_dict[col].tolist()
				for val in vals:
					## concatenate new column and set flag
					new_col = ''.join(['[', vcombo_name, ']', col, ':', control_value, '-', val])
					design_df = pd.concat([design_df, pd.Series([''] * design_df.shape[0], name=new_col)], axis=1)
					design_df.loc[temp_design_matrix.index[temp_design_matrix[col] == val], new_col] = '1'
	## sort dataframe
	if 'GENOTYPE' in contrast_cols:
	    design_df = design_df.sort_values('GENOTYPE' + contrast_cols)
	else:
	    design_df = design_df.sort_values('GENOTYPE' + contrast_cols)

	return design_df

def build_design_table(summary_df, cmp_cols, control=None):
	"""
	Automatically build design table that has basic comparison groups.
	"""
	print('... Building design table')
	## intialize design dataframe
	design_df = summary_df[['GENOTYPE', 'REPLICATE'] + cmp_cols].copy()
	#cmp_cols = ['GENOTYPE'] + cmp_cols
	## create a dict for col: vals
	col_dict = {}
	for col in cmp_cols:
		col_dict[col] = pd.unique(summary_df[col])
	## iterate thru columns to compare
	#for col in cmp_cols:
	for col in col_dict:
		if len(col_dict[col]) < 2: 
			## nothing to contrast with
			continue	
		## get value combination of other column descriptors
		other_cols = copy.copy(cmp_cols)
		other_cols.remove(col)
		n_other_cols = len(other_cols)
		other_vals = [col_dict[ocol] for ocol in other_cols]
		other_vals_combos = utils.makeListProduct(other_vals)
		## iterate thru with other column value fixed
		for vcombo in other_vals_combos:
			## get rows matching the column value combo
			df2 = summary_df[sum([summary_df[other_cols[k]] == vcombo[k] \
								for k in range(n_other_cols)]) == n_other_cols]
			vcombo_name = '-'.join([':'.join([str(other_cols[k]), str(vcombo[k])]) \
									for k in range(n_other_cols)])
			## assess diff column types
			if col == 'GENOTYPE' and control is not None:
				vals = col_dict[col].tolist()
				vals.remove(control)
				for val in vals:
					## concatenate new column and set flag
					new_col = ''.join(['[',vcombo_name,']',col,':',control,'-',val])
					design_df = pd.concat([design_df, pd.Series(['']*design_df.shape[0], name=new_col)], axis=1)
					design_df.loc[df2.index[df2[col] == control], new_col] = '0'
					design_df.loc[df2.index[df2[col] == val], new_col] = '1'
			else:
				vals = col_dict[col].tolist()
				for val in vals:
					## concatenate new column and set flag
					new_col = ''.join(['[', vcombo_name, ']', col, ':', '-', str(val)])
					design_df = pd.concat([design_df, pd.Series([''] * design_df.shape[0], name=new_col)], axis=1)
					design_df.loc[df2.index[df2[col] == control], new_col] = '0'
					design_df.loc[df2.index[df2[col] == val], new_col] = '1'
	## sort dataframe
	design_df = design_df.sort_values(cmp_cols)
	return design_df


def save_dataframe(filepath, df, df_cols=None):
	"""
	Save dataframe of quality assessment.
	"""
	if not filepath.endswith('.xlsx'):
		filepath += '.xlsx'
	df.to_excel(filepath, index=False, columns=df_cols)


if __name__ == '__main__':
	main(sys.argv)
