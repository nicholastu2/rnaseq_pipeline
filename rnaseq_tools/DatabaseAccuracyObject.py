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
            print('Cannot find .git/FETCH_HEAD in database_files. If this is a new, or newly cloned, directory, pull from the remote.')
        except AttributeError:
            print('.git/FETCH_HEAD is empty. Make a commit and try again.')
        # set accuracyCheckFilename (expecting to be overwritten by @property method below when needed)
        self.accuracy_check_output_file = self.setAccuracyCheckFilename()

    def setAccuracyCheckFilename(self):
        time_stamp = str(self.year_month_day) + '_' + utils.hourMinuteSecond()
        return os.path.join(self.reports, 'database_accuracy_check_' + time_stamp + '.txt')

    @staticmethod
    def subdirectoryReport(subdirectory_name, database_assessment_object, subdir_filepath_list):
        """

        """
        print('Checking %s column name formatting and entries' %subdirectory_name)
        specs_website = 'https://github.com/BrentLab/database_files/wiki'
        dba = database_assessment_object
        dba.setAccuracyCheckFilename()
        with open(dba.accuracy_check_output_file, 'a') as subdirectory_report:
            subdirectory_report.write('Checking %s for adherence to specifications found at: %s\n' % (subdirectory_name, specs_website))
            subdirectory_report.write('Last update (likely git pull) to directory: %s\n\n' % dba.last_git_change)

            for subdirectory_filepath in subdir_filepath_list:
                subdirectory_report.write('Checking %s:\n' % subdirectory_filepath)
                # extract dictionaries of inconsistencies in column names and rows
                col_inconsistencies_dict, row_inconsistencies_dict = dba.checkColumns(dba.specification_dict[subdirectory_name],
                                                                                      subdirectory_filepath, dba.logger)
                # check filename
                filename_check = dba.checkFileName(dba.specification_dict[subdirectory_name]['filename_regex'],
                                                   subdirectory_filepath)
                # check the format of the filename
                if not isinstance(filename_check, bool):
                    subdirectory_report.write('\tThe filename %s does not adhere to the specifications. Please correct.\n\n' % filename_check)
                # check column headings
                subdirectory_report.write('\tThe items below are column headings in a given sheet that do not match the specifications.\n')
                for spec_column, sheet_column in col_inconsistencies_dict.items():
                    subdirectory_report.write('\tThe specification is: %s, the sheet column is: %s\n' % (spec_column, sheet_column))
                subdirectory_report.write('\n\tThe items below are numbered by row (eg 0: inductionDelay means a problem in row 0 (first row) in column inductionDelay):\n')
                for row_index, column_heading in row_inconsistencies_dict.items():
                    subdirectory_report.write('\tRow %s has an inconsistency in column %s\n' % (row_index, column_heading))
                subdirectory_report.write('\n\n')
            subdirectory_report.write('\n\n')

    def fullReport(self):
        """
            The intent is for this to be used to generate a full report on the entire database. However, any number of
            subdirectories may be passed up to all of the subdirectories in database_files
        """
        for subdirectory_name, subdirectory_path_list in self.database_dict.items():
            self.subdirectoryReport(subdirectory_name, self, subdirectory_path_list)

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

    @staticmethod
    def checkFileName(filename_regex, subdirectory_filepath):
        """
            check that filename matches filename_regex.
            :param filename_regex: see DatabaseAccuracy Object filename_regexes
            :param subdirectory_filepath: filename of a individual sheet in a given database_files subdirectory
        """
        subdirectory_basename = os.path.basename(subdirectory_filepath)
        if not re.match(filename_regex, subdirectory_basename):  # TODO: ERROR CHECKING PATH
            return subdirectory_basename
        else:
            return True

    @staticmethod
    def checkColumns(subdirectory_specs_dict, subdirectory_filepath, logger = None):
        """
            check column heading names and entries in each row/column for adherence to the specs at:
            https://github.com/BrentLab/database_files/wiki
            :param subdirectory_specs_dict: see constructor. In the case of bioSample, you would pass db.specification_dict['bioSample']
            :param subdirectory_filepath: path to a sheet in a given subdirectory (eg a bioSample .xslx)
            :param logger: reference to a logger. Default is None
            :return: colname_inconsistencies_dict, a dict in structure {specification_heading: nearest_match_to_heading, ...}
                     row_inconsistencies_dict, a dict in structure {row_index: column_with_inconsistent_entry, ...}
        """
        # see :return: statement for structure
        colname_inconsistencies_dict = {}
        row_inconsistencies_dict = {}
        # list to store inappropriately formated column names
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
                        if logger:
                            logger.info('Column name not found in specs: %s' %(column_name))
                        nearest_match = difflib.get_close_matches(column_name, subdirectory_specs_dict['column_specs_dict'].keys())[0]
                        colname_inconsistencies_dict.setdefault(nearest_match, column_name)
                        print('\tCannot check %s in %s. Either the format of the column is incorrect, or it is not in the specifications_dictionary.\n'
                              '\tThe rest of this column could not be checked. Correct the column name, and re-run.' %(column_name, subdirectory_filepath))
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
        boolean_format = r"TRUE|FALSE|%r|%r|%i|%i" %(True, False, True, False)
        dna_alphabet_format = r"^[ACGT]+$"
        capital_lower_underscore_digit_format = r"^[a-zA-Z_\d]+$"
        capital_underscore_digit_format = r"[A-Z_\d]+"

        flood_media_options = r"None|PBS|SCGal|SCGlu"
        biosample_treatment = r"estradiol|mockEstradiol|conditionShift|glucoseTo2%|EtOH|Estradiol|PBS|" \
                              r"DMEM.30C.CO2.cAMP|RPMI.30C.CO2.cAMP|RPMI.37C.CO2.cAMP|YPD.30C.CO2.cAMP|" \
                              r"YPD.37C.CO2.cAMP|DMEM.37C.CO2|DMEM.30C.cAMP|DMEM.37C.cAMP|YPD.30C.cAMP|" \
                              r"YPD.37C.cAMP|RPMI.30C.cAMP|RPMI.37C.cAMP|RPMI.30C.CO2|RPMI.37C.CO2|YPD.30C.CO2|" \
                              r"YPD.37C.CO2"
        rna_prep_method = r"DirectZol|RiboPure0.5x|RiboPure0.25x|RiboPure0.125x|ComboA|ComboB|TRIzol"
        ribosomal_band_shape_options = r"straight|smile|NA"
        sequencer_model_options = r"NextSeq|MiSeq|MiniSeq"
        flowcell_options = r"V3|Standard|Nano|MiniSeq|HighOutput|MidOutput"
        purpose_options = r"Rebalancing|spikein|fullRNASeq|fullDNASeq|fullChIPSeq"
        s2cDNA_options = r"E7420L|SolexaPrep"

        bioSample_regex = r"^bioSample_[A-Z]\.[A-Z]+_\d+\.\d+\.\d+.[csvxlsx]+$"
        bioSample_column_dict = {'harvestDate': date_format,
                                 'harvester': name_format,
                                 'bioSampleNumber': int_format,
                                 'experimentDesign': capital_lower_underscore_digit_format,
                                 'experimentObservations': capital_lower_underscore_digit_format,
                                 'strain': capital_lower_underscore_digit_format,
                                 # TODO: LIST OF SPECIFIC POSSIBILITIES?
                                 'genotype': capital_underscore_digit_format,
                                 'floodmedia': flood_media_options,
                                 'inductionDelay': int_format,
                                 'treatment': biosample_treatment,
                                 'timePoint': int_format,
                                 'replicate': int_format}

        rnaSample_regex = r"^rnaSample_[A-Z]\.[A-Z]+_\d+\.\d+\.\d+.[csvxlsx]+$"
        rnaSample_column_dict = {'harvestDate': date_format,
                                 'harvester': name_format,
                                 'bioSampleNumber': int_format,
                                 'rnaDate': date_format,
                                 'rnaPreparer': name_format,
                                 'rnaSampleNumber': int_format,
                                 'rnaPrepMethod': rna_prep_method,
                                 # TODO: the next line in the specs is rnaPrepProtocol -- not in sheets
                                 'roboticRNAPrep': boolean_format,
                                 'RIBOSOMAL_BAND': boolean_format,
                                 'RIBOSOMAL_BAND_SHAPE': ribosomal_band_shape_options,
                                 'SMALL_RNA_BANDS': boolean_format,
                                 'RIN': one_through_ten_format}

        s1cDNASample_filename_regex = r"^s1cDNASample_[A-Z]\.[A-Z]+_\d+\.\d+\.\d+.[csvxlsx]+$"
        s1cDNASample_column_dict = {'rnaDate': date_format,
                                    'rnaPreparer': name_format,
                                    'rnaSampleNumber': int_format,
                                    's1cDNADate': date_format,
                                    's1cDNAPreparer': name_format,
                                    's1cDNASampleNumber': int_format,
                                    'polyAIsolationProtocol': r"%r" %"NEBNextPoly(A)E7490L",
                                    's1Protocol': r"E7420",
                                    'roboticS1Prep': boolean_format,
                                    's1PrimerSeq': r"^[ACGT]+$|random"}

        s2cDNASample_filename_regex = r"^s2cDNASample_[A-Z]\.[A-Z]+_\d+\.\d+\.\d+.[csvxlsx]+$"
        s2cDNASample_column_dict = {'s1cDNADate': date_format,
                                    's1cDNAPreparer': name_format,
                                    's1cDNASampleNumber': int_format,
                                    's2cDNADate': date_format,
                                    's2cDNAPreparer': name_format,
                                    's2cDNASampleNumber': int_format,
                                    's2cDNAProtocol': s2cDNA_options,
                                    'PooledSecondStrand': boolean_format,
                                    'roboticS2Prep': boolean_format}

        library_filename_regex = r"^library_[A-Z]\.[A-Z]+_\d+\.\d+\.\d+.[csvxlsx]+$"
        library_column_dict = {'s2cDNADate': date_format,
                               's2cDNAPreparer': name_format,
                               's2cDNASampleNumber': int_format,
                               'libraryDate': date_format,
                               'libraryPreparer': name_format,
                               'librarySampleNumber': int_format,
                               'index1Name': int_format,
                               'index1Sequence': dna_alphabet_format,
                               'index2Name': r"^SIC_Index_\d+$",
                               'index2Sequence': dna_alphabet_format,
                               'libraryProtocol': 'E7420',
                               'roboticLibraryPrep': boolean_format}

        fastqFilename_filename_regex = r"^fastqFiles_[A-Z]\.[A-Z]+_\d+\.\d+\.\d+.[csvxlsx]+$"
        fastqFilename_column_dict = {'libraryDate': date_format,
                                     'libraryPreparer': name_format,
                                     'librarySampleNumber': int_format,
                                     'runNumber': int_format,
                                     'laneNumber': int_format,
                                     'sequencerModel': sequencer_model_options,
                                     'flowcellType': flowcell_options,
                                     'purpose': purpose_options,
                                     'tapestationConc': float_format,
                                     'volumePooled': float_format,
                                     'readsObtained': int_format,
                                     'fastqFileName': capital_lower_underscore_digit_format}
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