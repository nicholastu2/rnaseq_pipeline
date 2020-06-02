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
        # set the database dictionary ({subdirectory: [list, of, files]} --> see DatabaseObject)
        self.setDatabaseDict()
        # if full_report passed in constructor
        try:
            if kwargs['full_report']:
                self.fullReport()  # Report output to report
        except KeyError:
            pass

        date_format = r"\d\d\.\d\d\.\d\d"
        name_format = r"[A-Z]\.[A-Z]*"
        sample_number_format = r"\d*"
        boolean_format = r"True|False"

        bioSample_regex = r"^bioSample_[A-Z]\.[A-Z]+_\d+\.\d+\.\d+.[csvxlsx]+$"
        bioSample_column_dict = {'harvestDate': date_format,
                                 'harvester': name_format,
                                 'bioSampleNumber': sample_number_format,
                                 'experimentDesign': r"[a-zA-Z_\d]*",
                                 'experimentObservations': r"[a-zA-Z_\d]*",
                                 'strain': r"[a-zA-Z_\d]*",  # TODO: LIST OF SPECIFIC POSSIBILITIES?
                                 'genotype': r"[A-Z_\d]*",
                                 'floodmedia': r"SCGal",
                                 'inductionDelay': r"\d*",
                                 'treatment': r"estradiol|mockEstradiol|conditionShift|glucoseTo2%|EtOH|Estradiol|PBS|DMEM.30C.CO2.cAMP|RPMI.30C.CO2.cAMP|RPMI.37C.CO2.cAMP|YPD.30C.CO2.cAMP|YPD.37C.CO2.cAMP|DMEM.37C.CO2",
                                 'timePoint': r"\d*",
                                 'replicate': r"\d*"}

        rnaSample_regex = r"^rnaSample_[A-Z]\.[A-Z]+_\d+\.\d+\.\d+.[csvxlsx]+$"
        rnaSample_column_dict = {'harvestDate': date_format,
                                 'harvester': name_format,
                                 'bioSampleNumber': sample_number_format,
                                 'rnaDate': date_format,
                                 'rnaPreparer': name_format,
                                 'rnaSampleNumber': sample_number_format,
                                 'rnaPrepMethod': r"DirectZol|RiboPure0.5x|RiboPure0.25x|RiboPure0.125x|ComboA|ComboB|TRIzol", #TODO: the next line in the specs is rnaPrepProtocol -- not in sheets
                                 'roboticRNAPrep': boolean_format,
                                 'RIBOSOMAL_BAND': boolean_format,
                                 'Ribosomal_Band_Shape': r"straight|smile|NA",
                                 'Small_RNA_Bands': boolean_format,
                                 'RIN': r"[1-9]|10"}

        s1cDNASample_filename_regex = r"^s1cDNASample_[A-Z]\.[A-Z]+_\d+\.\d+\.\d+.[csvxlsx]+$"
        s1cDNASample_column_dict = {'rnaDate': date_format,
                                    'rnaPreparer': name_format,
                                    'rnaSampleNumber': sample_number_format,
                                    's1cDNADate': date_format,
                                    's1cDNAPrepaper': name_format,
                                    's1cDNASampleNumber': sample_number_format,
                                    'PolyAIsolationProtocol': r"NEBNextPoly(A)E7490L",
                                    's1Protocol': r"E7420",
                                    'roboticS1Prep': r"True|False",
                                    's1PrimerSeq': r"[ACGT]+|random"}

        s2cDNASample_filename_regex = r"^s2cDNASample_[A-Z]\.[A-Z]+_\d+\.\d+\.\d+.[csvxlsx]+$"
        s2cDNASample_column_dict = {'s1cDNADate': date_format,
                                    's1cDNAPrepaper': name_format,
                                    's1cDNASampleNumber': sample_number_format,
                                    's2cDNADate': date_format,
                                    's2cDNAPrepaper': name_format,
                                    's2cDNASampleNumber': sample_number_format,
                                    's2cDNAProtocol': r"E7420|SolexaPrep",
                                    'PooledSecondStrand': r"True|False",
                                    'rotobics2Prep': r"True|False"}

        library_filename_regex = r"^library_[A-Z]\.[A-Z]+_\d+\.\d+\.\d+.[csvxlsx]+$"
        library_column_dict = {'s2cDNADate': date_format,
                               's2cDNAPrepaper': name_format,
                               's2cDNASampleNumber': sample_number_format,
                               'libraryDate': date_format,
                               'libraryPreparer': name_format,
                               'librarySampleNumber': sample_number_format,
                               'index1Name': r"\d*",
                               'index1Sequence': r"[ACGT]+",
                               'index2Name': r"SIC_Index_\d+",
                               'index2Sequence': r"[ACGT]+",
                               'libraryProtocol': 'E7420',
                               'roboticLibraryPrep': r"True|False"}

        fastqFilename_filename_regex = r"^fastqFiles_[A-Z]\.[A-Z]+_\d+\.\d+\.\d+.[csvxlsx]+$"
        fastqFilename_column_dict = {'libraryDate': date_format,
                                     'libraryPreparer': name_format,
                                     'librarySampleNumber': sample_number_format,
                                     'runNumber': sample_number_format,
                                     'laneNumber': sample_number_format,
                                     'sequencerModel': r"NextSeq|MiSeq|MiniSeq",
                                     'flowcellType': r"V3|Standard|Nano|MiniSeq|HighOutput|MidOutput",
                                     'purpose': r"Rebalancing|spikein|fullRNASeq|fullDNASeq|fullChIPSeq",
                                     'tapestationConc': r"\d+\.\d+|\d+",
                                     'volumePooled': r"\d+\.\d+|\d+",
                                     'readsObtained': r"\d+",
                                     'fastqFileName': r"[a-zA-Z_\d]*"}
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
        # set last_git_change
        self.last_git_change = self.getLastGitChange(self.database_files)
        # set accuracyCheckFilename (expecting to be overwritten by @property method below when needed)
        time_stamp = str(self.year_month_day) + '_' + utils.hourMinuteSecond()
        self.accuracy_check_output_file = os.path.join(self.reports, 'database_accuracy_check_' + time_stamp + '.txt')

    def setAccuracyCheckFilename(self):
        time_stamp = str(self.year_month_day) + '_' + utils.hourMinuteSecond()
        self.accuracy_check_output_file = os.path.join(self.reports, 'database_accuracy_check_' + time_stamp + '.txt')

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
                                                                                      subdirectory_filepath)
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

    @staticmethod
    def fullReport(database_dict, database_assessment_object):
        """
            The intent is for this to be used to generate a full report on the entire database. However, any number of
            subdirectories may be passed up to all of the subdirectories in database_files
        """
        dba = database_assessment_object
        dba.setAccuracyCheckFilename()
        for subdirectory_name, subdirectory_path_list in database_dict.items():
            database_assessment_object.subdirectoryReport(subdirectory_name, dba, subdirectory_path_list)

    @staticmethod
    def getLastGitChange(database_code_repo):
        """

        """
        git_head_file = os.path.join(database_code_repo, '.git/FETCH_HEAD')
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
    def checkColumns(subdirectory_specs_dict, subdirectory_filepath):
        """
            check column heading names and entries in each row/column for adherence to the specs at:
            https://github.com/BrentLab/database_files/wiki
            :param subdirectory_specs_dict: see constructor. In the case of bioSample, you would pass db.specification_dict['bioSample']
            :param subdirectory_filepath: path to a sheet in a given subdirectory (eg a bioSample .xslx)
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
                        nearest_match = difflib.get_close_matches(column_name, subdirectory_specs_dict['column_specs_dict'].keys())[0]
                        colname_inconsistencies_dict.setdefault(nearest_match, column_name)
                        print('\tCannot check %s in %s. Either the format of the column is incorrect, or it is not in the specifications_dictionary.\n'
                              '\tThe rest of this column could not be checked. Correct the column name, and re-run.' %(column_name, subdirectory_filepath))
                        skip_columns.append(column_name)
                else:
                    if not re.match(column_specs_regex, column_entry):
                        row_inconsistencies_dict.setdefault(str(index), []).append(column_name)

        return colname_inconsistencies_dict, row_inconsistencies_dict
