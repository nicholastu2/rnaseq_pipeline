from rnaseq_tools import utils
from rnaseq_tools import StandardData_tools
import pandas as pd
import os
import sys
import re
import time
import getpass  # see https://www.saltycrane.com/blog/2011/11/how-get-username-home-directory-and-hostname-python/


class StandardData:
    """
    parent class of rnaseq_pipeline. Creates rnaseq_pipeline directory in $USER if it (or any part) do not exist
    and stores all paths to resource directories and files.
       Does not check contents of all required directories for completeness (especially critical for genome_files), #TODO fix this so it does
       only that they exist.
       loads rnaseq_pipeline package level configuration (both main config as well as logger config in $USER/rnaseq_pipeline/config
           Recall that OrganismData, a child of StandardData, loads the .ini config files in each organism in $USER/rnaseq_pipeline/genome_files/
    """
    def __init__(self, expected_attributes=None, *args, **kwargs):
        """
        initialize StandardDataFormat with arbitrary number of keyword arguments.
        :param expected_attributes: a list of other attributes (intended use is for sub classes to add expected attributes)
        :param kwargs: arbitrary number/length keyword arguments. key = value will be set as class attributes
        """
        self.self_type = 'StandardData'
        # list of expected/standard attribute names
        self._attributes = ['lts_rnaseq_data', 'lts_sequence', 'lts_align_expr', 'scratch_database_files',
                            'scratch_sequence', 'opt_genome_files', 'user_rnaseq_pipeline', 'genome_files',
                            'reports', 'query', 'sbatch_log', 'log_dir', 'log_file', 'job_scripts', 'rnaseq_tmp',
                            'query_sheet_path', 'raw_count_path', 'log2_cpm_path', 'norm_counts_path', 'config_file',
                            'experiment_dir', 'email', 'fastq_path', 'strandness', 'run_number', 'output_dir',
                            'align_count_path', 'experiment_columns']

        self._run_numbers_with_zeros = {641: '0641', 647: '0647', 648: '0648', 659: '0659', 673: '0673', 674: '0674',
                                        684: '0684',
                                        731: '0731', 748: '0478', 759: '0759', 769: '0769', 773: '0773', 779: '0779'}

        # set year_month_day
        self.year_month_day = utils.yearMonthDay()
        # This is to extend _attributes if a class extends StandardData
        if isinstance(expected_attributes, list):
            self._attributes.extend(expected_attributes)
        # get user name and set as _user
        self._user = getpass.getuser()
        # set attributes entered by keyword on instantiation, warn user if keyword entered in instantiation not in _attributes
        kwargs['config_file'] = '/opt/apps/labs/mblab/software/rnaseq_pipeline/1.0/config/rnaseq_pipeline_config.ini'
        StandardData_tools.setAttributes(self, self._attributes, kwargs)
        # load config file
        utils.configure(self, self.config_file, self.self_type)
        # set interactive to false if not already set
        if not hasattr(self, 'interactive'):
            self.interactive = False           # TODO: NEED TO CHECK THIS -- CHILDREN MAY NOT OVERWITE!
        # create standard directory structure in /scratch/mblab/$USER (this path will be stored as self.scratch_rnaseq_pipeline)
        self.standardDirectoryStructure()

        # automatic actions to perform on certain attributes when instantiated
        if hasattr(self, 'query_sheet_path'):
            # set attribute 'query_df' to store the standardizedQuerySheet (and out to rnaseq_tmp)
            setattr(self, 'query_df', self.standardizeQuery(self.query_sheet_path))
            self.writeStandardizedQueryToTmp()  # TODO: sometime necessary to pass file btwn python and R. May be able to re-do this with package subprocess and pass stored variable rather than write?
        if hasattr(self, 'raw_counts_path'):    #               Alternatively, write to cluster tmp (give option to save)?
            # set attribute raw_counts_df to store raw_counts
            setattr(self, 'raw_counts_df',
                    pd.read_csv(self.raw_counts_path))  # TODO: decide if you actually need to do this

        # TODO: CREATE LOGGER
        # create instance of logger. This will be parent of each StandardData child logger and will read config file in rnaseq_pipeline/config
        #logger_configuration_file = ''  when logger config file is passed, enter here
        self.logger = utils.createLogger(self.log_file, __name__)  # this will create a logger with default settings without a config file -- see utils
        self.logger.debug("the year month day is %s" % self.year_month_day)

    def standardDirectoryStructure(self):
        """
        checks for and creates if necessary the expected directory structure in /scratch/mblab/$USER/rnaseq_pipeline
        """
        # first, create pipeline directory if dne
        setattr(self, 'user_scratch', os.path.join(self.mblab_scratch, self._user))
        setattr(self, 'user_rnaseq_pipeline', '{}/rnaseq_pipeline'.format(self.user_scratch))
        if not os.path.exists(self.user_rnaseq_pipeline):
            utils.mkdirp(self.user_rnaseq_pipeline)

        # check for the directories to be soft linked from /scratch/mblab/mblab.shared (self.mblab_shared). soft link if not, setattr in either case
        mblab_shared_dirs = ['scratch_sequence', 'database_files']
        self.softLinkAndSetAttr(self, mblab_shared_dirs, self.mblab_shared, self.user_rnaseq_pipeline)

        if not self.interactive:
            # check for directories to be soft linked from /lts/mblab/Crypto/rnaseq_pipeline (self.lts_rnaseq_data)
            lts_dirs_to_softlink = ['lts_align_expr', 'lts_sequence']
            self.softLinkAndSetAttr(self, lts_dirs_to_softlink, self.lts_rnaseq_data, self.user_rnaseq_pipeline)

        # unzip genome files from /lts/mblab/Crypto/rnaseq_data/1.0/genome_files to self.user_rnaseq_pipeline
        setattr(self, 'genome_files', os.path.join(self.user_rnaseq_pipeline, 'genome_files'))
        if not os.path.exists(self.genome_files):
            genome_files_full_path = os.path.join(self.lts_rnaseq_data, self.pipeline_version, 'genome_files.zip')
            cmd = 'unzip {} -d {}'.format(genome_files_full_path, self.user_rnaseq_pipeline)
            utils.executeSubProcess(cmd)

        # next, make directories if dne
        process_directories = ['reports', 'query', 'sbatch_log', 'log/%s' % self.year_month_day, 'job_scripts',
                               'rnaseq_tmp']
        for directory in process_directories:
            # store path
            path = os.path.join(self.user_rnaseq_pipeline, directory)
            # this will only create the path if it dne
            utils.mkdirp(path)
            # set attr to directory (the names in process_directories) unless log, which is treated specially
            if not directory == 'log/%s' % self.year_month_day:
                setattr(self, directory, path)
            else:
                # distinguish the log directory ($USER/rnaseq_pipeline/log)
                self.log_dir = 'log/%s' % self.year_month_day
                # from the daily log file ($USER/rnaseq_pipeline/log/<year-month-day>)
                self.log_file = os.path.join(self.log_dir, '%s.log' % self.year_month_day)

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
            # store the path (either it exists or it will after the if not block) to the directory in intended_dir_path)
            path = os.path.join(intended_dir_path, directory)
            # check if directory exists in the intended_dir_path
            if not os.path.exists(path):
                # if it does not, soft link ln -s origin_dir_path/directory to intended_dir_path/directory
                cmd = 'ln -s {}/{} {}'.format(origin_dir_path, directory, path)
                utils.executeSubProcess(cmd)
            # set attribute named directory (from for loop above) that points towards variable path which stores intended_dir_path/directory
            setattr(object_instance, directory, path)

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
        timestr = time.strftime("%Y%m%d_%H%M%S")
        standardized_query_name = '{}_standardized_query.csv'.format(timestr)
        standardized_query_path = os.path.join(self.rnaseq_tmp, standardized_query_name)
        self.query_df.to_csv(standardized_query_path, index=False)
        setattr(self, 'standardized_query_path', os.path.join(self.rnaseq_tmp, standardized_query_path))

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
        self.logger.debug('the extracted value, prior to running return_with_leading_zero is: %s' % extracted_value)
        return_with_leading_zero = check_leading_zero
        if return_with_leading_zero:  # TODO: casting this to an int is a bit ugly -- may be a point of weakness
            if int(extracted_value) in self._run_numbers_with_zeros:
                return str(self._run_numbers_with_zeros[int(extracted_value)])

        return extracted_value
