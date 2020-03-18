#!/usr/bin/env python
import os
import re
from abc import ABC, abstractmethod
import numpy as np
from itertools import combinations, product
import yaml
import sys
import pandas as pd
import getpass # see https://www.saltycrane.com/blog/2011/11/how-get-username-home-directory-and-hostname-python/
import time
#used in StandardDataFormat


def decompose_status2bit(n):
    """
    Decompose the bit status
    """
    if n == 0:
        return None
    decomp = []
    for i in np.arange(np.floor(np.log2(n)), -1, -1):
        if n - 2 ** i >= 0:
            decomp.append(i)
            n -= 2 ** i
    return decomp


def make_combinations(lst):
    """
    Make all possible replicate combinations
    """
    if len(lst) < 2:
        return [lst]
    combo = []
    for i in range(len(lst), 1, -1):
        for s in combinations(lst, i):
            combo.append(sorted(s))
    return combo


def makeListProduct(lst):
    """
    Make all possible value combination, with each value coming from a each list. Support up to 10 lists.
    """
    if len(lst) == 0:
        return None
    elif len(lst) == 1:
        return list(product(lst[0]))
    elif len(lst) == 2:
        return list(product(lst[0], lst[1]))
    elif len(lst) == 3:
        return list(product(lst[0], lst[1], lst[2]))
    elif len(lst) == 4:
        return list(product(lst[0], lst[1], lst[2], lst[3]))
    elif len(lst) == 5:
        return list(product(lst[0], lst[1], lst[2], lst[3], lst[4]))
    else:
        sys.exit('ERROR: The number of contrasts is greater than 10. This length beyond what I can handle')


def load_config(json_file):
    """
    Load configuration file (JSON) for QC thresholding and scoring
    """
    with open(json_file) as json_data:
        d = yaml.safe_load(json_data)
    return d


def parse_gtf(filename):
    """
    Convert gtf into dictionary for filling bed fields
    """
    bed_dict = {}
    reader = open(filename, 'r')
    for line in reader.readlines():
        if not line.startswith('#'):
            ## get info from each line
            line_split = line.split("\t")
            chrm, ltype, strand, annot = line_split[0], line_split[2], line_split[6], line_split[8]
            coords = [int(x) for x in line_split[3:5]]
            gene_id = re.findall(r'gene_id "(.+?)";', annot)[0]
            ## fill dictionary
            if gene_id not in bed_dict.keys():
                bed_dict[gene_id] = {'chrm': chrm, 'strand': strand, 'coords': [0, 0]}
            if ltype in 'start_codon':
                if strand == '+':
                    bed_dict[gene_id]['coords'][0] = min(coords)
                elif strand == '-':
                    bed_dict[gene_id]['coords'][1] = max(coords)
            elif ltype in 'stop_codon':
                if strand == '+':
                    bed_dict[gene_id]['coords'][1] = max(coords)
                elif strand == '-':
                    bed_dict[gene_id]['coords'][0] = min(coords)
    reader.close()
    return bed_dict


def parse_gff3(filename):
    """
    Convert gff3 into dictionary for filling bed fields
    """
    bed_dict = {}
    reader = open(filename, 'r')
    for line in reader.readlines():
        if not line.startswith('#'):
            ## get info from each line
            line_split = line.split("\t")
            chrm, ltype, strand, annot = line_split[0], line_split[2], line_split[6], line_split[8]
            coords = [int(x) for x in line_split[3:5]]
            if ltype in 'gene':
                gene_id = re.findall(r'Name=(.+?);', annot)[0]
                ## fill dictionary
                bed_dict[gene_id] = {'chrm': chrm, 'strand': strand, 'coords': [min(coords), max(coords)]}
    reader.close()
    return bed_dict


def mkdir_p(path_to_directory):
    """
    Function to create a directory. Equivalent to mkdir -p in bash. Will create directory if the path does not already exist.
    :param path_to_directory: path to a directory to be created if DNE already
    """
    if not os.path.exists(path_to_directory):
        try:
            os.makedirs(path_to_directory)
        except OSError:
            print("Directory {} cannot be created. See mkdir_p function call "
                  "in the script from which this error occurred".format(path_to_directory))



def check_dir(d):
    if not d.endswith('/'):
        d += '/'
    return d


def addForwardSlash(path):
    """
    add forward slash to path
    :param path: a path, typically to the lowest subdir prior to adding a file
    :returns: the path with a trailing slash added
    """
    if not path.endswith('/'):
        path = path + '/'
    return path


# the below functions strip file extensions, including in cases where there are two eg fq.gz
def fileBaseName(file_name):
    """
    strip everything after first . in path (may be too extreme -- make sure there are no 'necessary' . in filepaths. However, please try to avoid putting a . in a filepath other than the extension(s)
    :param file_name: the base filename of a file (not a path!)
    :returns: the filename stripped of all extensions
    """
    # https://stackoverflow.com/a/46811091
    if '.' in file_name:
        separator_index = file_name.index('.')
        base_name = file_name[:separator_index]
        return base_name
    else:
        return file_name


def pathBaseName(path):
    """
    strips directory structure and extensions. eg /path/to/file.fastq.gz --> file
    :param path: a filepath that includes a directory structure and file (see description)
    :returns: the basename of the file stripped of the directory structure and extensions
    """
    # This gets the basename of a given path, and then strips all file extensions (even if multiple). see fileBaseName
    file_name = os.path.basename(path)
    return fileBaseName(file_name)


def checkCSV(file):
    """
    test whether a given file is a .csv or something else
    :param file: a file that ends with .csv, .xlsx, .tsv, etc.
    :returns: True if .csv, false otherwise
    """
    # test whether a given file is a .csv or .xlsx
    if re.search('\.csv', file):
        return True
    else:
        return False


def getDirName(path):
    """
    get the directory name of the directory one level up from the end of the path provided. eg /path/to/file.fastq.gz returns 'to'
    :param path: a path to a directory or file
    :returns: the name of the directory one level up from the final level of path. see description for example
    """
    if os.path.split(path)[1] == "":
        dirname = os.path.dirname(path)
    else:
        dirname = os.path.basename(path)

    return dirname

class FileWriter:
    pass


class SlurmJobscriptWriter:
    pass

# TODO: input config of experiment as json in dictionary style, then take arguments from that
class CreateDesignMatrixColumns:
    """
        The intention of this is to create a ExperimentDesign class that is easily extensible to the different varieties of experiments
        we do in the lab. This is a rough (and explicit, as opposed to creating an interface and extending it) draft.
        Currently it is only set up for Zev timecourse experiments
    """

    def __init__(self, experimental_column_headers, contrast_column_header, quality_summary_df, control_value):
        # store cmd line input as class attributes
        self.experimental_column_headers = experimental_column_headers
        self.contrast_column_header = contrast_column_header
        self.quality_summary_df = quality_summary_df
        self.design_df_seed = self.createDesignMatrixSeed()  # remove columns other than experimental + contrast columns from sample_summary_df
        self.control_value = control_value
        # create experimental controls and contrast groups
        self.experimental_conditions_dict = {}
        self.experimental_conditions_iterable = self.createExperimentalConditionTuples()  # tuples in form ('GCN4', 'Estradiol') in case of [GENOTYPE, TREATMENT]
        self.contrast_conditions = self.createContrastConditions()  # this is without the control condition
        # Using the attributes, create the design matrix
        self.design_df = self.completeDesignMatrix()

    ### end constructor

    # TODO: include crypto/wt experiments. re-write as general 'experiment' class. extend for various types of experiments.

    def createDesignMatrixSeed(self):
        """
            parse out the experimental_condition_columns + contrast_columns
            of the quality_summary_df and return. This will be the beginning of the design matrix
            Returns:
        """
        columns_of_interest = self.experimental_column_headers + self.contrast_column_header

        if not 'GENOTYPE' in columns_of_interest:
            columns_of_interest.insert(0, 'GENOTYPE')
        if not 'REPLICATE' in columns_of_interest:
            columns_of_interest.insert(1, 'REPLICATE')
        if not 'FASTQFILENAME' in columns_of_interest:
            columns_of_interest.insert(2, 'FASTQFILENAME')

        return self.quality_summary_df[columns_of_interest]

    ### end createDesignMatrixSeed()

    def createExperimentalConditionTuples(self):
        """
            Returns: Tuples created by taking cartesian product of items in the unique values of the experimental_columns in the design_df
        """
        for column_heading in self.experimental_column_headers:
            self.experimental_conditions_dict.setdefault(column_heading, []).extend(
                list(pd.unique(self.design_df_seed[column_heading])))

        return product(*self.experimental_conditions_dict.values())

    ### end createExperimentalConditionsTuples()

    def createContrastConditions(self):
        """
            Returns: List of unique values from the contrast condition column, minus the control_value
        """
        contrast_conditions = list(np.unique(self.design_df_seed[self.contrast_column_header].values))
        contrast_conditions.remove(self.control_value)

        return contrast_conditions

    ### end createContrastConditions()

    def addDesignMatrixColumn(self, experimental_condition_tuple, contrast_value):
        """
            Creates a column heading in style explained in the FILES section of the rnaseq_pipeline github wiki
            Briefly, the heading will be [experimental conditions]contrast condition (eg [GENETYPE:GCN4 & TREATMENT:ESTRADIOL]TIMEPOINT:-1 & 10)
            Args:
            Returns:
        """
        # append to string while next() for some iterable to create column
        column_heading = ''
        # dictionaries to store column: value to use in identifying the rows in the design_df to add 0s and 1s to
        control_column_identifiers = {}
        contrast_column_identifiers = {}
        for i in range(len(self.experimental_column_headers)):
            # build new column header for the design table
            column_heading = column_heading + ' & ' + self.experimental_column_headers[i] + ':' + str(
                experimental_condition_tuple[i])
            # column: value to dictionaries
            contrast_column_identifiers.setdefault(self.experimental_column_headers[i], []).append(
                experimental_condition_tuple[i])
            control_column_identifiers.setdefault(self.experimental_column_headers[i], []).append(
                experimental_condition_tuple[i])
        # clean up column heading, add contrast group
        column_heading = '[' + column_heading[3:] + ']' + self.contrast_column_header[0] + ':' + str(
            self.control_value) + ' & ' + str(contrast_value)
        # add contrast column: value to dictionaries
        contrast_column_identifiers.setdefault(self.contrast_column_header[0], []).append(contrast_value)
        control_column_identifiers.setdefault(self.contrast_column_header[0], []).append(self.control_value)
        # create new column in df
        self.design_df_seed[column_heading] = ""
        # fill the column where appropriate with 0s and 1
        self.fillDesignMatrixColumn(column_heading, control_column_identifiers, contrast_column_identifiers)

    ### end createDesignMatrixColumnHeading()

    def fillDesignMatrixColumn(self, column_heading, control_column_identifiers, contrast_column_identifiers):
        """
            Fill the appropriate column with either a zero or 1
            Args:
            Returns:
        """
        # create boolean masks to identify which row is the control and which is the contrast. The dictionaries passed in the arguments are
        # used to filter, which returns 1s (boolean T) where there is a match. If the sum across the rows is equal to the number of keys in
        # the dictionary, then this is a positive match for a row in either the control or contrast group and will be assigned a 0/1 in the loop
        # below
        mask_control = self.design_df_seed.isin(control_column_identifiers).sum(axis=1) == len(
            control_column_identifiers)
        mask_contrast = self.design_df_seed.isin(contrast_column_identifiers).sum(axis=1) == len(
            contrast_column_identifiers)
        # fill row in design_df based on masks above
        for index, row in self.design_df_seed.iterrows():
            if mask_control[index] == True:
                self.design_df_seed.loc[index, column_heading] = '0'
            if mask_contrast[index] == True:
                self.design_df_seed.loc[index, column_heading] = '1'
            else:
                pass

    ### end fillDesignMatrixColumnHeading()

    def completeDesignMatrix(self):
        """
            Append design matrix columns to the design matrix seed
            Args:
            Returns:
        """
        for experimental_condition_tuple in self.experimental_conditions_iterable:
            for contrast_value in self.contrast_conditions:
                self.addDesignMatrixColumn(experimental_condition_tuple, contrast_value)
        # TODO: this needs to be cleaned up -- self.design_df should be modified in place while design_df_seed should be left alone.
        return self.design_df_seed

    ### end completeDesignMatrix()


### end CreateDesignMatrix

class StandardDataFormat:
    def __init__(self, **kwargs):
        """
        initialize StandardDataFormat with arbitrary number of keyword arguments.
        :param kwargs: arbitrary number/length keyword arguments. key = value will be set as class attributes
        """
        # list of expected attribute names
        self._attributes = ['query_sheet_path', 'raw_count_path', 'norm_counts_path', 'align_expr_path',
                            'sequence_path', 'database_files_path', 'genome_files', 'tmp_dir']
        self._user = getpass.getuser()
        self.tmp_dir = '/scratch/mblab/{}/rnaseq_tmp'.format(self._user)
        mkdir_p(self.tmp_dir) # make sure this is accessible to the class when reformat utils
        self.setAttributes(kwargs)
        if hasattr(self, 'query_sheet_path'):
            # set attribute 'query_df' to store the standardizedQuerySheet
            setattr(self, 'query_df', self.standardizeQuery(self.query_sheet_path))
            self.writeStandardizedQueryToTmp()
        if hasattr(self, 'raw_counts_path'):
            # set attribute raw_counts_df to store raw_counts
            setattr(self, 'raw_counts_df', pd.read_csv(self.raw_counts_path))

    def setAttributes(self, input_dict):
        for key, value in input_dict.items():
            if key not in self._attributes:
                print("{} not in expected attributes. This is not a problem unless you expect the data corresponding to "
                      "the values in {} to be automatically standardized")
            setattr(self, key, value)

    @staticmethod
    def standardizeQuery(df_path, prefix='', suffix='_read_count.tsv', only_sample_col=False):
        """
        convert a dataframe containing sample info to a 'standard form' -- capitalized column headings and FASTQFILENAME
        converted to SAMPLE with appropriate prefix and suffix replacing sequence/run_####_samples/ and .fastq.gz
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
            fastq_basename = pathBaseName(fastq_file_path)
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
                    prefix = addForwardSlash(prefix)
                df.loc[index, 'FASTQFILENAME'] = prefix + fastq_basename + suffix

        # rename fastqfilename column to 'sample'
        df.rename(columns={'FASTQFILENAME': 'COUNTFILENAME'}, inplace=True)

        # only_sample_col provides method of returning the filenames from FASTQFILENAME as a series
        if only_sample_col:
            return df['COUNTFILENAME']
        # else, return the entire dataframe
        else:
            return df

    def writeStandardizedQueryToTmp(self):
        """
        write self.query_df to self.tmp_dir. format is date_time_standardized_query.csv eg 20200312_115103_standardized_query.csv
        also sets attribute to this path
        """
        timestr = time.strftime("%Y%m%d_%H%M%S")
        standardized_query_name = '{}_standardized_query.csv'.format(timestr)
        standardized_query_path = os.path.join(self.tmp_dir, standardized_query_name)
        self.query_df.to_csv(standardized_query_path, index=False)
        setattr(self, 'standardized_query_path', os.path.join(self.tmp_dir, standardized_query_path))

    @staticmethod
    def countsPerMillion(raw_count_path, output_FULL_path):
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

def genotypeCheck(standard_data):
    """
    function to perform genotype check
    :param standard_data: an instance of StandardDataFormat with attributes query_sheet_path, raw_count_path, log2_cpm,
                          output_dir, experiment_columns
    :returns: None. However, outputs in subdir of align_counts a series of IGV screenshots and histograms corresponding to the genes in the
    """

    # verify standard_data has correct attributes and the paths to the query_sheet and raw_counts exit TODO: update these checks
    if not (hasattr(standard_data, 'standardized_query_path') and os.path.exists(standard_data.standardized_query_path)):
        sys.exit('the StandardDataFormat instance does not have attribute query_sheet_path')
    elif not (hasattr(standard_data, 'raw_count_path') and os.path.exists(standard_data.raw_count_path)):
        sys.exit('the StandardDataFormat instance does not have attribute raw_count_path')
    elif not (hasattr(standard_data, 'log2_cpm_path') and os.path.exists(standard_data.log2_cpm_path)):
        sys.exit('the StandardDataFormat instance does not have attribute log2_cpm_path')
    elif not (hasattr(standard_data, 'output_dir') and os.path.exists(standard_data.output_dir)):
        sys.exit('the StandardDataFormat instance does not have attribute output_dir')
    elif not standard_data.experiment_columns:
        print('no experiment columns provided. By default, genotype, timepoint, treatment will be used. Do you wish to continue (y/n)?\n')
        answer = input()
        if answer == 'n' | answer == 'No' | answer == N:
            sys.exit('relaunch script with appropriate experiment columns')
        setattr(standard_data, 'experiment_columns', ['genotype', 'timepoint', 'treatment'])

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


class QualityAssessmentObject(ABC):
    def __init__(self, source_dir, output_dir):
        """
        constructor
        :param source_dir: where to find the alignment and count files
        :param output_dir: the directory into which to deposit the quality_assessment
        :param columns: list of variables to use as column
        """
        self.source_dir = source_dir
        self.output_dir = output_dir
        super().__init__()

    @abstractmethod
    def createQualityAssessmentDataFrame(self):
        pass

    @abstractmethod
    def saveDataFrame(self):
        pass
