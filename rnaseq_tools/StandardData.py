from rnaseq_tools import utils
import pandas as pd
import os
import sys
import re
import time
import getpass # see https://www.saltycrane.com/blog/2011/11/how-get-username-home-directory-and-hostname-python/
import configparser

class StandardData:
    def __init__(self, expected_attributes = None, *args, **kwargs):
        """
        initialize StandardDataFormat with arbitrary number of keyword arguments.
        :param expected_attributes: a list of other attributes (intended use is for sub classes to add expected attributes)
        :param kwargs: arbitrary number/length keyword arguments. key = value will be set as class attributes
        """
        # list of expected/standard attribute names
        self._attributes = ['lts_rnaseq_data', 'lts_sequence', 'lts_align_expr', 'scratch_database_files',
                            'scratch_sequence', 'opt_genome_files', 'user_rnaseq_pipeline', 'align_expr',
                            'genome_files', 'reports', 'query', 'sbatch_log', 'log', 'job_scripts', 'rnaseq_tmp',
                            'query_sheet_path', 'raw_count_path', 'norm_counts_path']

        # This is to extend _attributes if a class extends StandardData
        if isinstance(expected_attributes, list):
            self._attributes.extend(expected_attributes)
        # load config file
        self.config_file = '/opt/apps/labs/mblab/software/rnaseq_pipeline/1.0/rnaseq_pipeline/config/rnaseq_pipeline_config.ini'
        self.configure()
        # get user name and set as _user
        self._user = getpass.getuser()
        # create standard directory structure in /scratch/mblab/$USER (this will be stored as self.scratch_rnaseq_pipeline)
        self.standardDirectoryStructure()
        # set attributes entered by keyword on instantiation, warn user if keyword entered in instantiation not in _attributes
        utils.setAttributes(self, self._attributes, kwargs)
        # automatic actions to perform on certain attributes when instantiated
        if hasattr(self, 'query_sheet_path'):
            # set attribute 'query_df' to store the standardizedQuerySheet (and out to rnaseq_tmp)
            setattr(self, 'query_df', self.standardizeQuery(self.query_sheet_path))
            self.writeStandardizedQueryToTmp()  # TODO: re-do this with package subprocess and store rather than write
        if hasattr(self, 'raw_counts_path'):
            # set attribute raw_counts_df to store raw_counts
            setattr(self, 'raw_counts_df', pd.read_csv(self.raw_counts_path))

    def configure(self):
        """
        reads and sets the attributes listed in ['StandardData'] in rnaseq_pipeline/config/rnaseq_pipeline_config.ini
        """
        # read config file
        config = configparser.ConfigParser()
        config.read(self.config_file)
        # set attributes for StandardData
        for key, value in config['StandardData'].items():
            setattr(self, key, value) # by default, values are read in as strings. Currently, all filepaths, so this is good

    def standardDirectoryStructure(self):
        """
        checks for and creates if necessary the expected directory structure in /scratch/mblab/$USER/rnaseq_pipeline
        """
        # first, create pipeline directory if dne
        setattr(self, 'scratch_rnaseq_pipeline', '/scratch/mblab/{}/rnaseq_pipeline'.format(self._user))
        if not os.path.exists(self.scratch_rnaseq_pipeline):
            utils.mkdirp(self.scratch_rnaseq_pipeline)

        # check for the directories to be soft linked from /scratch/mblab/mblab.shared (self.mblab_shared). soft link if not, setattr in either case
        mblab_shared_dirs = ['scratch_sequence', 'database_files']
        self.softLinkAndSetAttr(self, mblab_shared_dirs, self.mblab_shared, self.scratch_rnaseq_pipeline)
        # check for directories to be soft linked from /lts/mblab/Crypto/rnaseq_pipeline (self.lts_rnaseq_data)
        lts_dirs_to_softlink = ['lts_align_expr', 'lts_sequence']
        self.softLinkAndSetAttr(self, lts_dirs_to_softlink, self.lts_rnaseq_data, self.scratch_rnaseq_pipeline)
        # next, make directories if dne
        process_directories = ['reports', 'query', 'sbatch_log', 'log', 'job_scripts', 'rnaseq_tmp']
        for directory in process_directories:
            # store path
            path = os.path.join(self.scratch_rnaseq_pipeline, directory)
            # this will only create the path if it dne
            utils.mkdirp(path)
            # set the attribute
            setattr(self, directory, path)

    @staticmethod
    def softLinkAndSetAttr(object_instance, list_of_dirs, origin_dir_path, intended_dir_path):
        """
        creates soft links and stores the path as attributes of object_instance
        :param object_instance: an instance of a given object (StandardData in this case)
        :param list_of_dirs: a list of directories that exists in origin_dir_path
        :param origin_dir_path: path to where the directories exist (actual data)
        :param intended_dir_path: the place where you'd like the directories soft linked
        """
        # loop through directories in list_of_dirs
        for directory in list_of_dirs:
            # check if directory exists in the intended_dir_path
            if not os.path.exists(os.path.join(intended_dir_path, directory)):
                # if it does not, store the path to intended_dir_path/directory as path
                path = os.path.join(intended_dir_path, directory)
                # soft link ln -s origin_dir_path/directory to intended_dir_path/directory
                os.system('ln -s {}/{} {}'.format(origin_dir_path, directory, path))
            # set attribute named directory (from for loop above) that points towards intended_dir_path/directory
            setattr(object_instance, directory, path)

    @staticmethod
    def standardizeQuery(df_path, prefix='', suffix='_read_count.tsv', fastq_filename_rename = 'COUNTFILENAME', only_filename_col=False):
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