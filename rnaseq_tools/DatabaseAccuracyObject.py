"""
   A child of DatabaseObject to verify conformance of database with the specifications at
   https://github.com/BrentLab/database_files/wiki

   usage: query =
"""
import os
import re
import difflib
from rnaseq_tools import utils
import subprocess
from rnaseq_tools.DatabaseObject import DatabaseObject


# TODO: more error handling in functions
class DatabaseAccuracyObject(DatabaseObject):
    def __init__(self, expected_attributes=None, **kwargs):
        """
            constructor
            :param **kwargs: unspecified # keyword arguments. The keywords that are currently handled, if entered:
                                logger_path = path to the directory in which to deposit the logger

        """
        # Call StandardData (parent class) constructor
        self._add_expected_attributes = []
        super(DatabaseAccuracyObject, self).__init__(self._add_expected_attributes, **kwargs)
        self.self_type = 'DatabaseAccuracyObject'
        # set DatabaseAccuracyObjectLogger
        self.logger = utils.createStandardObjectChildLogger(self, __name__)
        # set the database dictionary ({subdirectory: [list, of, files]} --> see DatabaseObject)
        self.setDatabaseDict()
        # if full_report passed in constructor
        try:
            if kwargs['full_report']:
                self.fullReport()  # Report output to report
        except KeyError:
            pass
        # create specification dict -- see class metadataSpecificationObject below this class
        self.specification_dict = metadataSpecificationObject().specification_dict
        # set last_git_change
        try:
            self.last_git_change = self.getLastGitChange()
        except FileNotFoundError:
            print(
                'Cannot find .git/FETCH_HEAD in database_files. If this is a new, or newly cloned, directory, pull from the remote.')
        except AttributeError:
            print('.git/FETCH_HEAD is empty. Make a commit and try again.')
        # set accuracyCheckFilename (expecting to be overwritten by @property method below when needed)
        self.accuracy_check_output_file = self.accuracyCheckFilename()
        self.key_column_dict = {"fastqFiles": ['libraryDate', 'libraryPreparer', 'librarySampleNumber'],
                                "library": ['libraryDate', 'libraryPreparer', 'librarySampleNumber'],
                                "s2cDNASample": ['s2cDNADate', 's2cDNAPreparer', 's2cDNASampleNumber'],
                                "s1cDNASample": ['s1cDNADate', 's1cDNAPreparer', 's1cDNASampleNumber'],
                                "rnaSample": ['rnaDate', 'rnaPreparer', 'rnaSampleNumber'],
                                "bioSample": ['harvestDate', 'harvester', 'bioSampleNumber']}

    def accuracyCheckFilename(self, descriptor = None):
        """
            create an accuracy check filename in rnaseq_pipeline/reports. All filenames will have format database_accuracy_check_optionsDescriptor_yearMonthDay.txt
            :params descriptor: default, none. Enter a single word to insert into filename after check_ before _yearMonthDay
            :returns: filepath for accuracy check output
        """
        # create beginning of path
        time_stamp = str(self.year_month_day)
        accuracy_check_filename = os.path.join(self.reports,'database_accuracy_check_')

        # if descriptor is passed, add it here
        if descriptor:
            accuracy_check_filename =  accuracy_check_filename + descriptor + '_'

        # append the timestamp and return
        return accuracy_check_filename + time_stamp + '.txt'

    def subdirectoryReport(self, subdirectory_name, subdir_filepath_list, key_column_only_report=False):
        """

        """
        print('Checking %s column name formatting and entries' % subdirectory_name)
        specs_website = 'https://github.com/BrentLab/database_files/wiki'
        with open(self.accuracy_check_output_file, 'a') as subdirectory_report:
            subdirectory_report.write('Checking %s for adherence to specifications found at: %s\n' % (subdirectory_name, specs_website))
            subdirectory_report.write('Last update (likely git pull) to directory: %s\n\n' % self.last_git_change)

            for subdirectory_filepath in subdir_filepath_list:
                self.logger.debug('Checking %s:\n' % subdirectory_filepath)
                # extract dictionaries of inconsistencies in column names and rows
                col_inconsistencies_dict, row_inconsistencies_dict = self.checkColumns(self.specification_dict[subdirectory_name], subdirectory_filepath)
                # check the format of the filename
                if not self.checkFileName(subdirectory_name, subdirectory_filepath):
                    subdirectory_report.write('\tThe filename %s does not adhere to the specifications. Please correct.\n\n' % os.path.basename(subdirectory_filepath))
                # check column headings
                lines_to_write = ['In sheet %s:\n\tThe items below are column headings in a given sheet that do not match ' \
                                 'the specifications (key and non-key, this should be fixed when found).\n' %subdirectory_filepath]
                for spec_column, sheet_column in col_inconsistencies_dict.items():
                    lines_to_write.append('\tThe specification is: %s, the sheet column is: %s\n' % (spec_column, sheet_column))
                lines_to_write.append('\n\tThe items below are numbered by row (eg 1: inductionDelay means a problem in row 1 of inductionDelay). If shortReport, only key columns are checked:\n')
                for row_index, column_heading in row_inconsistencies_dict.items():
                    # if short_report flag == True, only write out if the column_heading is a key column
                    subdir_key_set = set(self.key_column_dict[utils.pathBaseName(utils.dirPath(subdirectory_filepath))])
                    current_column_heading_set = set(column_heading)
                    # determine if column heading is in key column set
                    key_set_diff_length = len(subdir_key_set - current_column_heading_set)
                    if not key_column_only_report or (len(subdir_key_set) != key_set_diff_length):
                        lines_to_write.append('\tRow %s has an inconsistency in column %s\n' % (row_index, column_heading))
                # if no columns found to have inconsistencies, remove the header line for this section from the lines_to_write list
                if lines_to_write[-1].endswith('only key columns are checked:\n'):
                    lines_to_write.pop(-1)
                # if no column headings are found to be inconsistent, don't write at all. otherwise, write out the lines
                if not lines_to_write[-1].endswith('this should be fixed when found).\n'):
                    lines_to_write.append('\n\n\n\n')
                    subdirectory_report.write(''.join(lines_to_write))

    def report(self, key_columns_only=False):
        """
            The intent is for this to be used to generate a full report on the entire database. However, any number of
            subdirectories may be passed up to all of the subdirectories in database_files
            :params key_columns_only: only check actual filename and key columns for adherence to specs
        """
        # remove old sheets from the same day if they exist
        if os.path.isfile(self.accuracy_check_output_file):
            remove_cmd = 'rm %s' %self.accuracy_check_output_file
            utils.executeSubProcess(remove_cmd)

        if key_columns_only:
            self.accuracy_check_output_file = self.accuracyCheckFilename('keyColumn')

        for subdirectory_name, subdirectory_path_list in self.database_dict.items():
            self.subdirectoryReport(subdirectory_name, subdirectory_path_list, key_columns_only)

    def getLastGitChange(self):
        """

        """
        git_head_file = os.path.join(self.database_files, '.git/FETCH_HEAD')
        if not os.path.isfile(git_head_file):
            raise FileNotFoundError('GitHeadNotFound')
        stat_output = subprocess.run('stat -c %y {}'.format(git_head_file), stdout=subprocess.PIPE, shell=True)
        try:
            last_modified_date = stat_output.stdout.strip()
        except AttributeError:
            print('stat_output from %s failed')
        else:
            if last_modified_date == '':
                raise AttributeError('NoStatOutput')
            else:
                return stat_output.stdout.strip()

    def checkFileName(self, subdirectory_name, subdirectory_filepath):
        """
            check that filename matches filename_regex.
            :param subdirectory_name: name of a subdirectory in database files
            :param subdirectory_filepath: filename of a individual sheet in a given database_files subdirectory
        """
        # error check filepath
        try:
            if not os.path.isfile(subdirectory_filepath):
                raise FileNotFoundError('NotValidPath - %s' %subdirectory_filepath)
        except FileNotFoundError:
            error_msg = '%s not a valid path' % subdirectory_filepath
            print(error_msg)
            self.logger.critical(error_msg)

        # test filename for adherence to filename specs
        filename_regex = self.specification_dict[subdirectory_name]['filename_regex']
        subdirectory_basename = os.path.basename(subdirectory_filepath)
        file_check_flag = True
        if not re.match(filename_regex, subdirectory_basename):  # TODO: ERROR CHECKING PATH
            file_check_flag = False

        return file_check_flag

    def checkColumns(self, subdirectory_specs_dict, subdirectory_filepath):
        """
            check column heading names and entries in each row/column for adherence to the specs at:
            https://github.com/BrentLab/database_files/wiki
            :param subdirectory_specs_dict: see constructor. In the case of bioSample, you would pass db.specification_dict['bioSample']
            :param subdirectory_filepath: path to a sheet in a given subdirectory (eg a bioSample .xslx)
            :param logger: reference to a logger. Default is None
            :return: colname_inconsistencies_dict, a dict in structure {specification_heading: nearest_match_to_heading, ...}
                     row_inconsistencies_dict, a dict in structure {row_index: column_with_inconsistent_entry, ...}
        """
        self.logger.info('path to sheet is %s' % subdirectory_filepath)
        # see :return: statement for structure
        colname_inconsistencies_dict = {}
        row_inconsistencies_dict = {}
        # list to store inappropriately formatted column names
        skip_columns = []

        # read in subdirectory_filepath as dataframe
        subdirectory_df = utils.readInDataframe(subdirectory_filepath)
        # loop over rows in dataframe
        for index, row in subdirectory_df.iterrows():
            # convert row into a dictionary {column: value, ...}
            row_dict = dict(row)
            for column_name, column_entry in row_dict.items():
                column_entry = str(column_entry)
                try:
                    column_specs_regex = subdirectory_specs_dict['column_specs_dict'][column_name]
                except KeyError:
                    if column_name not in skip_columns:
                        if self.logger:
                            self.logger.info('Column name not found in specs: %s' % column_name)
                            self.logger.info('row for offending column is: %s' % row)
                        nearest_match = difflib.get_close_matches(column_name, subdirectory_specs_dict['column_specs_dict'].keys())[0]
                        colname_inconsistencies_dict.setdefault(nearest_match, column_name)
                        print('\tCannot check %s in %s. Either the format of the column is incorrect, or it is not in the specifications_dictionary.\n'
                              '\tThe rest of this column could not be checked. Correct the column name, and re-run.' % (column_name, subdirectory_filepath))
                        skip_columns.append(column_name)
                else:
                    if not re.match(column_specs_regex, column_entry):
                        row_inconsistencies_dict.setdefault(str(index), []).append(column_name)

        return colname_inconsistencies_dict, row_inconsistencies_dict


class metadataSpecificationObject:
    def __init__(self):
        """
            A dictionary to store the specifications of the Brentlab rnaseq metadata database. specs published at
            https://github.com/BrentLab/database_files/wiki
        """
        date_format = r"^\d\d\.\d\d\.\d\d$"
        name_format = r"^[A-Z]\.[A-Z]+$"
        int_format = r"^[ 0-9-]+$"
        one_through_ten_format = r"[1-9]|10"
        float_format = r"^\d+\.\d+$|^\d+$"
        # capitalized string, python True/False, and 0/1
        boolean_format = r"TRUE|FALSE|%r|%r|%i|%i" % (True, False, True, False)
        dna_alphabet_format = r"^[ACGT]+$"
        atleast_one_non_space_format = r"^(?!\s*$).+"
        capital_underscore_digit_format = r"[A-Z_\d]+"

        flood_media_options = r"None|PBS|SCGal|SCGlu"
        biosample_medium = r"DMEM|YPD|RPMI|PBS|mouseSerum"
        marker_format = r"NAT|G418|NaN|nan"
        rna_prep_method = r"DirectZol|RiboPure0.5X|RiboPure0.25X|RiboPure0.125X|ComboA|ComboB|TRIzol"
        ribosomal_band_shape_options = r"straight|smile|NA"
        sequencer_model_options = r"NextSeq|MiSeq|MiniSeq"
        flowcell_options = r"V3|Standard|Nano|MiniSeq|HighOutput|MidOutput"
        purpose_options = r"Rebalancing|spikein|fullRNASeq|fullDNASeq|fullChIPSeq"
        s2cDNA_options = r"SolexaPrep|E7420L\d+.\d+X|E7420L"

        bioSample_regex = r"^bioSample_[A-Z]\.[A-Z]+_\d+\.\d+\.\d+.[csvxlsx]+$|bioSample_\d+"
        bioSample_column_dict = {'harvestDate': date_format or r"Retrofit_\d+",
                                 'harvester': name_format or r"Retrofit_\d+",
                                 'bioSampleNumber': int_format or r"Retrofit_\d+",
                                 'experimentDesign': atleast_one_non_space_format or r"Retrofit_\d+",
                                 'experimentObservations': atleast_one_non_space_format or r"Retrofit_\d+",
                                 'strain': atleast_one_non_space_format or r"Retrofit_\d+",
                                 'genotype1': capital_underscore_digit_format,
                                 'genotype2': capital_underscore_digit_format or None,
                                 'perturbation1': capital_underscore_digit_format or r"Retrofit_\d+",
                                 'perturbation2': capital_underscore_digit_format or r"Retrofit_\d+",
                                 'floodmedia': flood_media_options or r"Retrofit_\d+",
                                 'inductionDelay': int_format or r"Retrofit_\d+",
                                 'medium': biosample_medium or None,
                                 'temperature': float_format or None,
                                 'atmosphere': "CO2" or None,
                                 'treatment': 'cAMP' or None,
                                 'other_conditions': 'pH7' or None,
                                 'timePoint': int_format or r"Retrofit_\d+",
                                 'replicate': int_format or r"Retrofit_\d+",
                                 'marker_1': marker_format or None,
                                 'marker_2': marker_format or None}

        rnaSample_regex = r"^rnaSample_[A-Z]\.[A-Z]+_\d+\.\d+\.\d+.[csvxlsx]+$|rnaSample_\d+"
        rnaSample_column_dict = {'harvestDate': date_format or r"Retrofit_\d+",
                                 'harvester': name_format or r"Retrofit_\d+",
                                 'bioSampleNumber': int_format or r"Retrofit_\d+",
                                 'rnaDate': date_format or r"Retrofit_\d+",
                                 'rnaPreparer': name_format or r"Retrofit_\d+",
                                 'rnaSampleNumber': int_format or r"Retrofit_\d+",
                                 'rnaPrepMethod': rna_prep_method or r"Retrofit_\d+",
                                 # TODO: the next line in the specs is rnaPrepProtocol -- not in sheets
                                 'roboticRNAPrep': boolean_format or r"Retrofit_\d+",
                                 'RIBOSOMAL_BAND': boolean_format or r"Retrofit_\d+",
                                 'RIBOSOMAL_BAND_SHAPE': ribosomal_band_shape_options or r"Retrofit_\d+",
                                 'SMALL_RNA_BANDS': boolean_format or r"Retrofit_\d+",
                                 'RIN': one_through_ten_format or r"Retrofit_\d+"}

        s1cDNASample_filename_regex = r"^s1cDNASample_[A-Z]\.[A-Z]+_\d+\.\d+\.\d+.[csvxlsx]+$|s1cDNASample_\d+"
        s1cDNASample_column_dict = {'rnaDate': date_format or r"Retrofit_\d+",
                                    'rnaPreparer': name_format or r"Retrofit_\d+",
                                    'rnaSampleNumber': int_format or r"Retrofit_\d+",
                                    's1cDNADate': date_format or r"Retrofit_\d+",
                                    's1cDNAPreparer': name_format or r"Retrofit_\d+",
                                    's1cDNASampleNumber': int_format or r"Retrofit_\d+",
                                    'polyAIsolationProtocol': r"E7490L_\d+\.\d+X|E7490L|catcher|None|Retrofit_\d+",
                                    's1cDNAProtocol': r"E7420L_\d+\.\d+X|E7420L|Retrofit_\d+|SuperScriptIII",
                                    'roboticS1Prep': boolean_format or r"Retrofit_\d+",
                                    's1PrimerSeq': r"^[ACGT]+$|random|Retrofit_\d+"}

        s2cDNASample_filename_regex = r"^s2cDNASample_[A-Z]\.[A-Z]+_\d+\.\d+\.\d+.[csvxlsx]+$|s2cDNASample_\d+"
        s2cDNASample_column_dict = {'s1cDNADate': date_format or r"Retrofit_\d+",
                                    's1cDNAPreparer': name_format or r"Retrofit_\d+",
                                    's1cDNASampleNumber': int_format or r"Retrofit_\d+",
                                    's2cDNADate': date_format or r"Retrofit_\d+",
                                    's2cDNAPreparer': name_format or r"Retrofit_\d+",
                                    's2cDNASampleNumber': int_format or r"Retrofit_\d+",
                                    's2cDNAProtocol': s2cDNA_options or r"Retrofit_\d+",
                                    'PooledSecondStrand': boolean_format or r"Retrofit_\d+",
                                    'roboticS2Prep': boolean_format or r"Retrofit_\d+"}

        library_filename_regex = r"^library_[A-Z]\.[A-Z]+_\d+\.\d+\.\d+.[csvxlsx]+$|library_\d+"
        library_column_dict = {'s2cDNADate': date_format or r"Retrofit_\d+",
                               's2cDNAPreparer': name_format or r"Retrofit_\d+",
                               's2cDNASampleNumber': int_format or r"Retrofit_\d+",
                               'libraryDate': date_format or r"Retrofit_\d+",
                               'libraryPreparer': name_format or r"Retrofit_\d+",
                               'librarySampleNumber': int_format or r"Retrofit_\d+",
                               'index1Name': r"Index\d_\d+",
                               'index1Sequence': dna_alphabet_format or r"Retrofit_\d+",
                               'index2Name': r"^SIC_Index2_\d+$|universal|Universal|Retrofit_\d+",
                               'index2Sequence': dna_alphabet_format or r"Retrofit_\d+",
                               'libraryProtocol': r"E7420L_\d+\.\d+X|E7420L" or r"Retrofit_\d+",
                               'roboticLibraryPrep': boolean_format or r"Retrofit_\d+"}

        fastqFilename_filename_regex = r"^fastqFiles_\d+\.\d+\.\d+.[csvxlsx]+$|fastq_\d+"
        fastqFilename_column_dict = {'libraryDate': date_format or r"Retrofit_\d+",
                                     'libraryPreparer': name_format or r"Retrofit_\d+",
                                     'librarySampleNumber': int_format or r"Retrofit_\d+",
                                     'runNumber': int_format or r"Retrofit_\d+",
                                     'laneNumber': int_format or r"Retrofit_\d+",
                                     'sequencerModel': sequencer_model_options or r"Retrofit_\d+",
                                     'flowcellType': flowcell_options or r"Retrofit_\d+",
                                     'purpose': purpose_options or r"Retrofit_\d+",
                                     'tapestationConc': float_format or None or r"Retrofit_\d+",
                                     'volumePooled': float_format or None or r"Retrofit_\d+",
                                     'readsObtained': int_format or None or r"Retrofit_\d+",
                                     'fastqFileName': atleast_one_non_space_format,
                                     'manualAudit': r"^(0|1)" or None,
                                     'manualStatus': r"[0-9]*" or None,
                                     'NOTES': r'[a-zA-Z]*' or None}

        # collect above specs into dictionary
        self.specification_dict = {
            'bioSample': {'filename_regex': bioSample_regex, 'column_specs_dict': bioSample_column_dict},
            'rnaSample': {'filename_regex': rnaSample_regex, 'column_specs_dict': rnaSample_column_dict},
            's1cDNASample': {'filename_regex': s1cDNASample_filename_regex,
                             'column_specs_dict': s1cDNASample_column_dict},
            's2cDNASample': {'filename_regex': s2cDNASample_filename_regex,
                             'column_specs_dict': s2cDNASample_column_dict},
            'library': {'filename_regex': library_filename_regex, 'column_specs_dict': library_column_dict},
            'fastqFiles': {'filename_regex': fastqFilename_filename_regex,
                           'column_specs_dict': fastqFilename_column_dict}}
