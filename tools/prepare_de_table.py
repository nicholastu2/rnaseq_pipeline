#!/usr/bin/env python
import sys
import argparse
import os.path
import pandas as pd
from rnaseq_tools import utils
from rnaseq_tools.ExperimentObject import CreateDesignMatrixColumns

def main(argv):
    parsed = parse_args(argv)
    if not os.path.exists(parsed.sample_summary):
        sys.exit('ERROR: %s does not exist.' % parsed.sample_summary)

    ## get conditions
    if isinstance(parsed.experimental_conditions, str):
        experimental_conditions = [c.strip() for c in parsed.experimental_conditions.split(' ')]
    else:
        experimental_conditions = parsed.experimental_conditions
    contrast_condition = parsed.contrast_conditions

    ## load sample summary
    summary_df = pd.read_excel(parsed.sample_summary)
    ## sift data based on group and quality
    summary_df = summary_df[(summary_df['MANUAL_AUDIT']==0)]
    ## prepare design table
    design_df = buildDesignTable(summary_df, experimental_conditions, contrast_condition, parsed.control_value)

    experiment_name = utils.fileBaseName(os.path.basename(parsed.sample_summary))
    experiment_name = experiment_name.replace('_quality_summary_2', '')
    design_table_filepath = os.path.join(parsed.output_dir, experiment_name + '_design_table')
    save_dataframe(design_table_filepath, design_df)

def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sample_summary', required=True,
                        help='[Required] <your_experiment>_quality_summary_2.xlsx -- the output of quality_assess_2.py')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='[Required] Where to deposit the design table for DE anlaysis. Suggested usage: the experiment directory. The name of the file will be the basename of your sample_summary + _design_table.xlsx>')
    parser.add_argument('--experimental_conditions', nargs='+', default='TREATMENT TIMEPOINT',
                        help='[Default TREATMENT TIMEPOINT] Experimental conditions that describe the sample are used to identify subgroups within each '
                             'genotype. No delimiter necessary if multiple conditions used.')
    parser.add_argument('-c', '--contrast_conditions', nargs='+', default='GENOTYPE',
                        help='Experimental conditions to describe the sample. List without delimiter if multiple eg GENOTYPE INDUCTIONDELAY. Default is currently GENOTYPE (crypto)')
    parser.add_argument('-cv',  '--control_value', default='CNAG_00000',
                        help='Control, e.g. CNAG_00000 for crypto or -1 for timecourse. Default is currently CNAG_00000')
    return parser.parse_args(argv[1:])

def buildDesignTable(summary_df, experimental_conditions, contrast_condition, control):
    """
    Automatically build design table that has basic comparison groups.
    """
    # see CreateDesignMatrixColumns class in ExperimentObject.py in rnaseq_tools
    experiment_design = CreateDesignMatrixColumns(experimental_conditions, contrast_condition, summary_df, control)
    return experiment_design.design_df

def save_dataframe(filepath, df, df_cols=None):
    """
    Save dataframe of quality assessment.
    """
    if not filepath.endswith('.xlsx'):
        filepath += '.xlsx'
    df.to_excel(filepath, index=False, columns=df_cols)


if __name__ == '__main__':
    main(sys.argv)
