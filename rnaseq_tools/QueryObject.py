from rnaseq_tools import utils
import pandas as pd
import os
import re
import sys
from rnaseq_tools import utils


class QueryObject:
    def __init__(self):
        pass

    @staticmethod
    def standardizeQuery(df_path, prefix='', suffix='_read_count.tsv', fastq_filename_rename='COUNTFILENAME',
                         only_filename_col=False):
        """
        convert a dataframe containing sample info to a 'standard form' -- capitalized column headings and FASTQFILENAME
        converted to COUNTFILENAME with appropriate prefix and suffix replacing sequence/run_####_samples/ and .fastq.gz
        :param df_path: path to a dataframe containing (at least) fastqFileName column
        :param prefix: a prefix to attach to fastqFileName column (eg '/lts/mblab/Crypto/rnaseq_data/align_expr') default none
        :param suffix: What to append to fastqFileName after stripping .fastq.gz. default _read_count.tsv
        :param fastq_filename_rename: name to give the column fastqFileName. default COUNTFILENAME
        :param only_filename_col: boolean, returns only the filename column
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
                if prefix:  # if a prefix other than align_expr/run_{} is passed, test if forward slash is present
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
        write self.query_df (standardized query_sheet) to self.rnaseq_tmp. format is date_time_standardized_query.csv
        eg 20200312_115103_standardized_query.csv also sets attribute to this path
        """
        timestr = ('%s_%s' %(utils.yearMonthDay(), utils.hourMinuteSecond()))
        standardized_query_name = '{}_standardized_query.csv'.format(timestr)
        standardized_query_path = os.path.join(self.rnaseq_tmp, standardized_query_name)
        self.query_df.to_csv(standardized_query_path, index=False)
        setattr(self, 'standardized_query_path', os.path.join(self.rnaseq_tmp, standardized_query_path))

    def extractValueFromStandardizedQuery(self, filter_column, filter_value, extract_column, check_leading_zero=False):
        """
        extract a value from a row (selected by filter_value) of self.query_df
        :param filter_column:
        :param filter_value:
        :param extract_column:
        :param check_leading_zero: if true, return the run number as a string with a leading 0 if it is in self._run_numbers_with_zeros
        :returns: a value extracted from a certain column of a certain row
        """
        row = self.query_df[self.query_df[filter_column] == filter_value]

        extracted_value = row[extract_column].values[0]
        if extracted_value == '0478':
            self.logger.debug('the extracted value, prior to running return_with_leading_zero is: %s' % extracted_value)
            self.logger.debug('the row is: %s' % row)
            self.logger.debug('the filter_value is: %s' % filter_value)
        return_with_leading_zero = check_leading_zero
        if return_with_leading_zero:  # TODO: casting this to an int is a bit ugly -- may be a point of weakness
            if int(extracted_value) in self._run_numbers_with_zeros:
                return str(self._run_numbers_with_zeros[int(extracted_value)])

        return extracted_value

    # TODO: TAKE CARE OF THIS IN QUERYDB OBJECT (CHILD OF STANDARDDATA POSSIBLY)
    # # automatic actions to perform on certain attributes when instantiated
    # if hasattr(self, 'query_sheet_path'):
    #     # set attribute 'query_df' to store the standardizedQuerySheet (and out to rnaseq_tmp)
    #     setattr(self, 'query_df', self.standardizeQuery(self.query_sheet_path))
    #     self.writeStandardizedQueryToTmp()  # TODO: figure out how to handle this -- either delete when exit, or don't write until necessary.