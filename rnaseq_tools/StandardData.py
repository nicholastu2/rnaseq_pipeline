from rnaseq_tools import utils
import pandas as pd
import os
import sys
import re
import time
import getpass # see https://www.saltycrane.com/blog/2011/11/how-get-username-home-directory-and-hostname-python/

class StandardData:
    def __init__(self, expected_attributes = None, **kwargs):
        """
        initialize StandardDataFormat with arbitrary number of keyword arguments.
        :param expected_attributes: a list of other attributes (intended use is for sub classes to add expected attributes)
        :param kwargs: arbitrary number/length keyword arguments. key = value will be set as class attributes
        """
        # list of expected/standard attribute names
        self._attributes = ['query_sheet_path', 'raw_count_path', 'norm_counts_path', 'align_expr_path',
                            'sequence_path', 'database_files_path', 'genome_files', 'tmp_dir']
        if isinstance(expected_attributes, list):
            self._attributes.extend(expected_attributes)
        # get $USER username
        self._user = getpass.getuser()
        # create rnaseq_tmp in $USER scratch space if it does not already exist and store the path
        self.tmp_dir = '/scratch/mblab/{}/rnaseq_tmp'.format(self._user)
        utils.mkdirp(self.tmp_dir) # make sure this is accessible to the class when reformat utils
        # set attributes entered by keyword on instantiation, warn user if keyword entered in instantiation not in _attributes
        utils.setAttributes(self, self._attributes, kwargs)
        # automatic actions to perform on certain attributes when instantiated
        if hasattr(self, 'query_sheet_path'):
            # set attribute 'query_df' to store the standardizedQuerySheet
            setattr(self, 'query_df', self.standardizeQuery(self.query_sheet_path))
            self.writeStandardizedQueryToTmp()
        if hasattr(self, 'raw_counts_path'):
            # set attribute raw_counts_df to store raw_counts
            setattr(self, 'raw_counts_df', pd.read_csv(self.raw_counts_path))

    @staticmethod
    def standardizeQuery(df_path, prefix='', suffix='_read_count.tsv', fastq_filename_rename = 'COUNTFILENAME', only_filename_col=False):
        """
        convert a dataframe containing sample info to a 'standard form' -- capitalized column headings and FASTQFILENAME
        converted to COUNTFILENAME with appropriate prefix and suffix replacing sequence/run_####_samples/ and .fastq.gz
        :param df_path: path to a dataframe containing (at least) fastqFileName column
        :returns: the dataframe with column variables cast to uppercase and fastqFileName converted to SAMPLE
        """
        df = pd.read_csv(df_path)
        # convert column headings to upper case
        df.columns = df.columns.str.upper()

        # regex to extract run_number, if needed
        regex = r"(?<=sequence\/run_)\d*"

        # loop through rows
        for index, row in df.iterrows():
            # replace fastqfilename one by one so as to extract the run number appropriately
            fastq_file_path = df.loc[index, 'FASTQFILENAME']
            fastq_basename = utils.pathBaseName(fastq_file_path)
            if prefix == 'align_expr/run_{}/':
                try:
                    run_number = re.search(regex, fastq_file_path)[0]
                except TypeError:
                    sys.exit(
                        'No run number found in the path provided in the \'FASTQFILENAME\' column. See utils function'
                        'standardizeSampleDataFrame')
                df.loc[index, 'FASTQFILENAME'] = prefix.format(run_number) + fastq_basename + suffix
            else:
                if prefix: # if a prefix other than align_expr/run_{} is passed, test if forward slash is present
                    prefix = utils.addForwardSlash(prefix)
                df.loc[index, 'FASTQFILENAME'] = prefix + fastq_basename + suffix

        # rename fastqfilename column to 'sample'
        df.rename(columns={'FASTQFILENAME': fastq_filename_rename}, inplace=True)

        # only_sample_col provides method of returning the filenames from FASTQFILENAME as a series
        if only_filename_col:
            return df[fastq_filename_rename]
        # else, return the entire dataframe
        else:
            return df

    def writeStandardizedQueryToTmp(self):
        """
        write self.query_df (standardized query_sheet) to self.tmp_dir. format is date_time_standardized_query.csv
        eg 20200312_115103_standardized_query.csv also sets attribute to this path
        """
        timestr = time.strftime("%Y%m%d_%H%M%S")
        standardized_query_name = '{}_standardized_query.csv'.format(timestr)
        standardized_query_path = os.path.join(self.tmp_dir, standardized_query_name)
        self.query_df.to_csv(standardized_query_path, index=False)
        setattr(self, 'standardized_query_path', os.path.join(self.tmp_dir, standardized_query_path))

    @staticmethod
    def countsPerMillion(raw_count_path, output_FULL_path):
        """
        submit raw_count_path to log2_cpm.R (in tools/)
        :param raw_count_path: path to output of raw_counts.py
        :param output_FULL_path: the full path (including the file and extension) of the output of log2_cpm.R. eg <experiment_name>_log2_cpm.csv
        """
        os.system('log2_cpm.R -r raw_count_path -o output_FULL_path')

    @staticmethod
    def userInputCorrectPath(message, object, attribute):
        # method to prompt user to correct path
        print(message)
        new_attr_value = input
        return setattr(object, attribute, new_attr_value)

    @staticmethod
    def userInputCorrectAttributeName(message, object, old_attribute_name):
        """
        usage: object = userInputCorrectAttributeName(...)
        """
        print(message + " ")
        new_attribute_name = input()
        setattr(object, new_attribute_name, getattr(object, old_attribute_name))
        delattr(object, old_attribute_name)
        # to use this, user will need to object = userInputCorrectAttributeName(...)
        return object