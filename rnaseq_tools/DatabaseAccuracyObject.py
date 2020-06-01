"""
   A child of DatabaseObject to verify conformance of database with the specifications at
   https://github.com/BrentLab/database_files/wiki

   usage: query =
"""
import pandas as pd
import os
import re
import sys
from rnaseq_tools import utils
import subprocess
from rnaseq_tools.DatabaseObject import DatabaseObject


# TODO: more error handling in functions
class DatabaseAccuracyObject(DatabaseObject):
    def __init__(self, **kwargs):
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
                self.fullReport() # Report output to report
        except KeyError:
            pass

    def bioSampleSpecifications(self):
        """

        """
        filename_regex= 'bioSample_J.PLAGGENBERG_05.30.19.csv'
        column_dict = {'harvestDate': None,
                       'harvester': None,
                       'bioSampleNumber': None,
                       'experimentDesign': None,
                       'experimentObservation': None,
                       'strain': None,
                       'genotype': None,
                       'floodmedia': None,
                       'inductionDate': None,
                       'treatment': None,
                       'timePoint': None,
                       'replicate': None}

    def rnaSampleSpecifications(self):
        """

        """
        filename_regex = 'rnaSample_J.PLAGGENBERG_05.30.19.csv'
        column_dict = {'harvestDate': None,
                       'harvester': None,
                       'bioSampleNumber': None,
                       'rnaDate': None,
                       'rnaPreparer': None,
                       'rnaSampleNumber': None,
                       'rnaPrepMethod': None,
                       'rnaPrepProtocol': None,
                       'roboticRNAPrep': None,
                       'Ribosomal_Band': None,
                       'Ribosomal_Band_Shape': None,
                       'Small_RNA_Bands': None,
                       'RIN': None}

    def s1cDNASampleSpecifications(self):
        """

        """
        filename_regex = 's1cDNASample_J.PLAGGENBERG_05.30.19.csv'
        column_dict = {'rnaDate': None,
                       'rnaPreparer': None,
                       'rnaSampleNumber': None,
                       's1cDNADate': None,
                       's1cDNAPrepaper': None,
                       's1cDNASampleNumber': None,
                       'PolyAIsolationProtocol': None,
                       's1Protocol': None,
                       'roboticS1Prep': None,
                       's1PrimerSeq': None}

    def s2cDNASampleSpecifications(self):
        """

        """
        filename_regex = 's2CDNASample_J.PLAGGENBERG_05.30.19.csv'
        column_dict = {'s1cDNADate': None,
                       's1cDNAPrepaper': None,
                       's1cDNASampleNumber': None,
                       's2cDNADate': None,
                       's2cDNAPrepaper': None,
                       's2cDNASampleNumber': None,
                       's2cDNAProtocol': None,
                       'PooledSecondStrand': None,
                       'rotobics2Prep': None}

    def librarySpecifications(self):
        """

        """
        filename_regex = 'library_J.PLAGGENBERG_05.30.19.csv'
        column_dict = {'s2cDNADate': None,
                       's2cDNAPrepaper': None,
                       's2cDNASampleNumber': None,
                       'libraryDate': None,
                       'libraryPreparer': None,
                       'librarySampleNumber': None,
                       'index1Name': None,
                       'index1Sequence': None,
                       'index2Name': None,
                       'index2Sequence': None,
                       'libraryProtocol': None,
                       'roboticLibraryPrep': None}

    def fastqFilesSpecifications(self):
        """

        """
        filename_regex = 'fastqFiles_J.PLAGGENBERG_05.30.19.csv'
        column_dict = {'libraryDate': None,
                       'libraryPreparer': None,
                       'librarySampleNumber': None,
                       'runNumber': None,
                       'laneNumber': None,
                       'sequencerModel': None,
                       'flowcellType': None,
                       'purpose': None,
                       'tapestationConc': None,
                       'volumnPooled': None,
                       'readsObtained': None,
                       'fastqFileName': None}

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
    def checkFileName(filename_regex, subdirectory_filename):
        """

        """
        raise NotImplementedError

    @staticmethod
    def checkColumns(column_specifications_dict, subdirectory_file_df):
        """

        """
        raise NotImplementedError

    @staticmethod
    def subdirectoryReport(subdirectory_specifications):
        """

        """
        # include "fields not included"
        raise NotImplementedError

    @staticmethod
    def fullReport(database_directory, subdirectory_list):
        """

        """
        raise NotImplementedError