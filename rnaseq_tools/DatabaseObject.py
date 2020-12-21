"""
   A class to interact with the brentlab rnaseq_database

   usage: query = DatabaseObject()
                  # keyword arguments of interest:
                       config_file = '/path/to/StandardDataObject_config_file'
                       interactive = True/False
                       And any number of database specific keywords, such as:
                           filter_json_path
                           database_df
                           etc (see constructor)
"""
import pandas as pd
import os
from rnaseq_tools import utils
from rnaseq_tools.StandardDataObject import StandardData

# TODO: more error handling in functions
class DatabaseObject(StandardData):
    def __init__(self, expected_attributes=None, **kwargs):
        """
            constructor
            :param **kwargs: unspecified # keyword arguments. The keywords that are currently handled, if entered:
                                logger_path = path to the directory in which to deposit the logger

            PLEASE NOTE: order of database_subdirectories is important. see comment in constructor

        """
        # Call StandardData (parent class) constructor
        self._add_expected_attributes = ['database_subdirectories', 'filter_json_path', 'database_df']
        super(DatabaseObject, self).__init__(self._add_expected_attributes, **kwargs)
        self.self_type = 'DatabaseObject'
        # set DatabaseObject logger
        self.logger = utils.createStandardObjectChildLogger(self, __name__)
        try:
            self.database_directory = kwargs['database_files']
        except KeyError:
            self.database_directory = self.database_files

        # set default database subdirectories. PLEASE NOTE: order is important here -- list in the order you wish them to merge in
        try:
            self.database_subdirectories = kwargs['database_subdirectories']
        except KeyError:
            #self.database_subdirectories = ['fastqFiles', 'library', 's2cDNASample', 's1cDNASample', 'rnaSample', 'bioSample'] # concat in reverse order
            self.database_subdirectories = ['bioSample', 'rnaSample', 's1cDNASample', 's2cDNASample', 'library', 'fastqFiles']

        # see setter setFilterJson()
        self.filter_json = None
        try:
            self.filter_json_path = kwargs['filter_json_path']
        except KeyError:
            self.filter_json_path = None

        # see setter setDatabaseDataframe()
        try:
            self.database_df = kwargs['database_df']
        except KeyError:
            self.database_df = None

        # see filterDatabaseDataframe()
        try:
            self.filtered_database_df = kwargs['filtered_database_df']
        except KeyError:
            self.filtered_database_df = None

        # data_dir_dict will store {database_subdirectory: [list, of, filepaths, in, each, subdir], ... }. See self.setDatabaseDict()
        self.database_dict = {}
        # see setter setDatabaseDict()
        self.concat_database_dict = {}
        # see setter setKeyColumns()
        self.database_key_columns = []

    def setDatabaseDict(self):
        """
            create dictionary of filepaths to the various types of metadata sheets
            :set_attribute data_dir_dict: structure {'subdirectory': [list, of, filenames]}
        """
        # associate each key (relevant subdirectories of database) with a list of of files (absolute path) in each key directory
        for subdirectory in self.database_subdirectories:
            # create a path database_files/subdirectory
            subdirectory_path = os.path.join(self.database_directory, subdirectory)
            # extract list of files in subdirectory_path (not recursive -- will return subdirs, but not their contents)
            subdirectory_files = utils.extractTopmostFiles(subdirectory_path)

            no_tmp_subdirectory_files = []
            # remove temporary files
            for file in subdirectory_files:
                basename = os.path.basename(file)

                if basename.startswith('~') or basename.startswith('._') or \
                        basename.startswith('.~') or basename.startswith('~$'):
                    pass
                else:
                    no_tmp_subdirectory_files.append(file)

            # test whether any of the key subdirectories of datadir are empty, throw error if so
            # TODO: make this error handling not so silly -- create error or handle more directly
            try:
                if len(no_tmp_subdirectory_files) == 0:
                    raise ValueError('NoSubdirFiles')
            except ValueError:  # TODO: include this in tests
                print("No files found in %s. "
                      "These files are necessary to creating the sample_summary." % subdirectory_path)
            else:
                # associate subdirectory (key) with subdirectory_files (value) in data_dir_dict
                self.database_dict.setdefault(subdirectory, []).extend(no_tmp_subdirectory_files)

    def setDatabaseDataframe(self, accuracy_check=False):
        """
            create joined data frame from the concatenated files in the subdirectories of the database_directory
            :param accuracy_check: boolean flag to indicate whether the purpose of concatenating the database is checking the string format accuracy. If true, name keys are not cast to uppper
        """
        # check that database_dict, concat_database_dict and database_key_columns exist
        if len(self.database_dict) == 0:
            self.setDatabaseDict()
        if len(self.concat_database_dict) == 0:
            self.setConcatDatabaseDict()
        if len(self.database_key_columns) == 0:
            self.setKeyColumns()
        # TODO: CAST ALL KEY COLUMNS TO UPPERCASE PRIOR TO MERGE
        # merge the first two (fastqFiles and Library) sets of data
        left_sheet = self.concat_database_dict[self.database_subdirectories[0]]
        right_sheet = self.concat_database_dict[self.database_subdirectories[1]]
        # if not an accuracy check (default) cast the name column (the second item in the key list) to upper case
        if not accuracy_check:
            left_sheet[self.database_key_columns[0][1]] = left_sheet[self.database_key_columns[0][1]].str.upper()
            right_sheet[self.database_key_columns[0][1]] = right_sheet[self.database_key_columns[0][1]].str.upper()
        # merge
        self.database_df = pd.merge(left_sheet, right_sheet, how='left', on=list(self.database_key_columns[0]))
        # merge the subsequent sheets on the columns identified in key_cols
        for i in range(1, len(self.database_subdirectories) - 1):
            # store the following for readability
            subdirectory = self.database_subdirectories[i + 1]
            database_key_column = self.database_key_columns[i]
            next_df = self.concat_database_dict[subdirectory]
            if not accuracy_check:
                next_df[database_key_column[1]] = next_df[database_key_column[1]].str.upper() # cast name column to upper
            # keep merging the next sheet to self.database_df
            self.database_df = pd.merge(self.database_df, self.concat_database_dict[subdirectory], how='left',
                                        on=list(database_key_column))

    def setConcatDatabaseDict(self):
        """
            creates concatenated dataframe from all files in a given list of paths to database subdirectories
            structure {subdirectory: concatenated_table_of_all_files_in_database/subdirectory/*}
        """
        # create dataframe from first file in file_list in given key(subdirectory) of database_dict
        for subdirectory, file_list in self.database_dict.items():
            self.concat_database_dict[subdirectory] = utils.readInDataframe(file_list[0])
            column_list = self.concat_database_dict[subdirectory].columns
            column_list = [column_header.strip() for column_header in column_list]
            self.concat_database_dict[subdirectory].columns = column_list
            # keep appending (cbind) dataframes to the bottom
            for file in file_list[1:]:
                # read in next file in list as next_sheet
                next_sheet = utils.readInDataframe(file)
                self.logger.debug('columns of %s are %s' %(file, next_sheet.columns))
                self.concat_database_dict[subdirectory] = self.concat_database_dict[subdirectory].append(next_sheet)
            # reset index so it is sequential
            self.concat_database_dict[subdirectory].reset_index(inplace=True, drop=True)

    def setKeyColumns(self):
        """
            create list of shared columns between successive pairs of keys in datadir_keys i.e. the columns which are
            shared between the sheets in fastqFiles and Library. See diagram of database at https://github.com/BrentLab/database_files/wiki

            PLEASE NOTE: Dictionaries are ordered in python 3.6+. key: directories need to be added to self.concat_database_dict in order you wish to merge
        """
        num_keys = len(self.database_subdirectories)
        # compare concat_data_base columns from one sheet to the next
        for i in range(num_keys - 1):
            # store the ith and ith+1 concat_database sheets
            left_database_columns = self.concat_database_dict[self.database_subdirectories[i]].keys()
            right_database_columns = self.concat_database_dict[self.database_subdirectories[i + 1]].keys()
            # store the key columns in list database_key_columns
            self.database_key_columns.append(left_database_columns.intersection(right_database_columns))

    def setFilterJson(self):
        """ TODO: error handling not working for invalid json. fix this.
            read in the filter json. See note in filterDatabaseDataframe()
            :raises: ValueError("NoFilterJson")
        """
        if self.filter_json_path is None:
            raise FileNotFoundError("NoFilterJson")
        else:
            try:
                self.filter_json = pd.read_json(self.filter_json_path, typ='series', dtype=False)
            except FileNotFoundError:
                print('problem with json')
                self.logger.error('No json present in filterDatabaseDataframe, or the formatting isn\'t recognized. Check the json -- try double quotes if using single -- and try again.')
            except ValueError:
                print('problem with json')
                self.logger.error('No json present in filterDatabaseDataframe, or the formatting isn\'t recognized. Check the json -- try double quotes if using single -- and try again.')

    def dropRowsIfEmptyFastqFilename(self):
        """
            drops rows if the fastqFilename entry is empty
        """
        try:
            self.concat_database_dict['fastqFiles'].dropna(subset=['fastqFileName'], inplace=True)
        except KeyError:
            self.logger.error('Unable to drop rows from concat_database_dict[\'fastqFiles\']. '
                              'Check that it exists as both a subdirectory and a dataframe')

    def filterDatabaseDataframe(self):
        """
            filter database_dataframe by json.
            Please note: this has limited capability. More complicated filtering should be done using a DatabaseObject
            and the pandas sql-like filtering commands
            see https://pandas.pydata.org/docs/getting_started/comparison/comparison_with_sql.html
        """
        try:
            if self.database_df is None:
                raise ValueError('DatabaseDataframeNotSet')
        except ValueError:
            self.setDatabaseDataframe()
        finally:
            if self.filter_json is None:
                self.setFilterJson()
                self.logger.debug(self.filter_json)
            # begin a string to store the query formula
            filter_str = '('
            # loop through columns in json query (i.e. 'timePoint' and 'treatment')
            for key, value in self.filter_json.items():
                filter_str = filter_str + '{} == {}'.format(key, value) + ' & '
            # once out of both loops, split the string on the final &, retain the substring before the split, and eliminate whitespace
            filter_str = filter_str.rsplit('&', 1)[0].strip() + ')'
            self.logger.debug('the filter created from the json_dict is %s' % filter_str)
            # use the filter_str formula to filter the dataframe
            self.filtered_database_df = self.database_df.query(filter_str)

    @staticmethod  #TODO: THIS NEEDS TO BE CHANGED SO THAT THE STANDARD FASTQFILENAME --> SAMPLE WITH JUST THE BASENAME OF THE SAMPLE, NO PATH OR EXTENSION
    def standardizeDatabaseDataframe(rnaseq_metadata_df, **kwargs):
        """
            convert a dataframe containing sample info to a 'standard form' -- capitalized column headings and FASTQFILENAME
            is just the sample name -- no path, no extension
            :param rnaseq_metadata_df: pandas dataframe of the rnaseq_metadata
            :param kwargs: arbitrary keyword arguments. provided to pass logger
            :returns: the dataframe with column variables cast to uppercase and fastqFileName converted to SAMPLE
        """
        try:
            # convert column headings to upper case
            rnaseq_metadata_df.columns = rnaseq_metadata_df.columns.str.upper()
        except AttributeError:
            print('standardizeDatabaseDataframe takes a dataframe, not a filepath, as an argument')

        # loop through rows
        for index, row in rnaseq_metadata_df.iterrows():
            # replace fastqfilename one by one so as to extract the run number appropriately
            fastq_file_path = rnaseq_metadata_df.loc[index, 'FASTQFILENAME']
            fastq_basename = utils.pathBaseName(fastq_file_path)
            rnaseq_metadata_df.loc[index, 'FASTQFILENAME'] = utils.pathBaseName(fastq_basename)

        return rnaseq_metadata_df

    @staticmethod
    def writeDatabaseDataframeToTmp(rnaseq_metadata_df, **kwargs):
        """
            write dataframe to file. Default write to current working directory with name rnaseq_metadata_date_time.csv
            eg rnaseq_metadata_20200312_115103.csv
            :param rnaseq_metadata_df: a dataframe -- expecting a database dataframe
            :param kwargs: arbitrary number of keyword arguments. Currently set to handle
                               path_to_directory: directory in which to write a file
                               filename: filename to replace rnaseq_metadata above. date/time will still be added, as will .csv
                                         DO NOT ENTER EITHER DATE/TIME OR .csv. eg filename = <my_experiment>_standardized
                                         would produce <my_experiment>_standardized_20200312_115103.csv
                                logger: a logger, likely the one created by DatabaseObject, but could pass any
            :returns: path to dataframe
        """
        timestr = ('%s_%s' % (utils.yearMonthDay(), utils.hourMinuteSecond()))
        # store path to directory. current working directory is default
        if kwargs['path_to_directory']:  # TODO error handling == this wont work
            path_to_directory = kwargs['path_to_directory']
        else:
            path_to_directory = os.getcwd()
        # store filename. rnaseq_metadata_yearmonthday_hourminutesecond.csv is default
        if kwargs['filename']:
            filename = kwargs['filename'] + '_%s' % timestr
        else:
            filename = 'rnaseq_metadata_%s.csv' % timestr
        # create path to dataframe
        dataframe_path = os.path.join(path_to_directory, filename)
        if kwargs['logger']:
            kwargs['logger'].debug('The dataframe path is: %s' % dataframe_path)
        # write dataframe
        rnaseq_metadata_df.to_csv(dataframe_path, index=False)
        # return path
        return dataframe_path

    @staticmethod
    def extractValueFromStandardizedQuery(rnaseq_metadata_df, filter_column, filter_value, extract_column, **kwargs):
        """
            extract a value from a row (selected by filter_value) of rnaseq_metadata_df. Currently set to handle only unique extractions (not lists)
            :param rnaseq_metadata_df: a rnaseq_metadata dataframe
            :param filter_column: a column on which to filter based on filter_value
            :param filter_value: a value in the filter_column on which to select a row
            :param extract_column: column from which to actually extract value (may not be same as filter_column)
            :param kwargs: arbitrary number of keyword arguments. Currently handles:
                               logger: a logger, probably based from DatabaseObject but could be any logger
                               leading_zero_dict: see StandardData._run_numbers_with_zeros
            :returns: a value extracted from a certain column of a certain row
        """
        # extract row of interest based on filter value
        row = rnaseq_metadata_df[rnaseq_metadata_df[filter_column] == filter_value]
        # extract value of interest from row
        extracted_value = row[extract_column].values[0]
        # if leading_zero_dict (see StandardData) passed to kwargs, if the extracted value is a runNumber in the leading_zero_dict, then return with a leading 0
        try:
            if kwargs['leading_zero_dict']:  # TODO: casting this to an int is a bit ugly -- may be a point of weakness
                if int(extracted_value) in kwargs['leading_zero_dict']:
                    return str(kwargs['leading_zero_dict'][int(extracted_value)])
        except KeyError:
            pass

        return extracted_value

    @staticmethod
    def uniqueKeys(key_columns, database_subdirectory_concat_df):
        """
            test whether the key columns in a given sheet are unique. This is an interactive method -- user enters values to stdin
            :param key_columns: key columns (see DatabaseObject.setKeyColumns())
            :param database_subdirectory_concat_df: see DatabaseObject.setConcatDatabaseDict()
        """
        # make a tuple of the key columns and store them as pandas series
        key_tuples = database_subdirectory_concat_df[key_columns].apply(tuple, axis=1)
        num_keys = key_tuples.size
        num_unique_keys = key_tuples.unique().size
        print(
            "\nThe number of unique keys is {}. The number of rows is {}. If these are equal, the keys are unique.".format(
                num_keys, num_unique_keys))
        if not num_keys == num_unique_keys:
            non_unique_rows = database_subdirectory_concat_df[database_subdirectory_concat_df[key_columns].duplicated()]
            print("\nThe following indicies are not unique:\n\t%s" % non_unique_rows)
            print("\nDo you want to continue? Enter y or n: ")
            user_response = input()
            if user_response == 'n':
                quit()

    # TODO: write coersion functions to coerce values in columns, correct known small typos (case, name mispells/misformats, etc)
    # def coerceAllCols(sheet):
    #     # converts columns specific in functions below to floats and datetime
    #     sheet = floatColCoerce(sheet)
    #
    #     sheet = datetimeColCoerce(sheet)
    #
    #     # forcing all to lowercase is not working in the search
    #     #sheet = strColCoerce(sheet)
    #
    #     return sheet
    #
    # def colCoerce(cols, sheet, dtype):
    #     # general function to coerce column to certain data type
    #     # Args: a list of columns to coerce (can be a subset of columns of a data frame. must be a list), a dataframe,
    #     #       and a datatype (see pandas .astype documentation for which datatypes are possible. didn't work for datetime, for example)
    #     # Return: the dataframe with the specified columns coerced
    #     coerce_cols = cols
    #
    #     for col in sheet.columns[sheet.columns.isin(coerce_cols)]:
    #         sheet[col] = sheet[col].astype(dtype)
    #     return sheet
    #
    # def floatColCoerce(sheet):
    #     # coerce columns specified below to int
    #     int_columns = ['librarySampleNumber', 'runNumber', 'tapestationConc', 'readsObtained', 'rnaSampleNumber',
    #                    'inductionDelay', 'replicate', 'timePoint']
    #
    #     sheet = colCoerce(int_columns, sheet, 'float')
    #
    #     return sheet
    #
    # def datetimeColCoerce(sheet):
    #     # coerce columns specified below to datetime
    #     datetime_columns = ['libraryDate', 's2cDNADate', 's1cDNADate', 'rnaDate', 'harvestDate']
    #
    #     for col in datetime_columns:
    #         sheet[col] = pd.to_datetime(sheet[col])
    #
    #     return sheet
    #
    # def strColCoerce(sheet):
    #     # coerce any column NOT in int_col or datetime_col to a string and make it lower case
    #
    #     int_columns = ['librarySampleNumber', 'tapestationConc', 'readsObtained', 'rnaSampleNumber','inductionDelay', 'replicate', 'timePoint']
    #     datetime_columns = ['libraryDate', 's2cDNADate', 's1cDNADate', 'rnaDate', 'harvestDate']
    #
    #     typed_cols = int_columns + datetime_columns
    #
    #     for col in sheet.columns:
    #         if col not in typed_cols:
    #             sheet[col] = sheet[col].astype(str).apply(lambda x: x.lower())
    #
    # return sheet
