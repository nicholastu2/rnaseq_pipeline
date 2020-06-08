import unittest
import os
from rnaseq_tools.DatabaseAccuracyObject import DatabaseAccuracyObject

class MyTestCase(unittest.TestCase):
    def test_setLastChange(self):
        """

        """
        test_database_files_path = '/home/chase/code/brentlab/rnaseq_pipeline/tests/test_data/database_files'
        function_last_change = DatabaseAccuracyObject.getLastGitChange(test_database_files_path)
        known_last_change = b'2020-05-30 08:56:23.755705986 -0500'

        self.assertEqual(known_last_change, function_last_change)

    def test_checkFileName(self):
        """

        """
        bioSample_path_correct = 'bioSample/bioSample_J.PLAGGENBERG_01.06.20.xlsx'
        biosample_basename_correct = os.path.basename(bioSample_path_correct)
        bioSample_regex = "^bioSample_[A-Z]\.[A-Z]*_\d\d\.\d\d\.\d\d\.[csvxlsx]*$"
        test_correct = DatabaseAccuracyObject.checkFileName(bioSample_regex, biosample_basename_correct)

        self.assertTrue(test_correct)

        bioSample_path_incorrect = 'bioSample/bioSample_jplag_01.06.20.xlsx'
        biosample_basename_incorrect = os.path.basename(bioSample_path_incorrect)

        test_incorrect = DatabaseAccuracyObject.checkFileName(bioSample_regex, biosample_basename_incorrect)

        self.assertEqual(test_incorrect, biosample_basename_incorrect)

    def test_checkColumns(self):
        db = DatabaseAccuracyObject(config_file='/home/chase/code/brentlab/rnaseq_pipeline/config/test_rnaseq_pipeline_config.ini')
        func_col_dict, func_row_dict = db.checkColumns(db.specification_dict['bioSample'], '/home/chase/code/brentlab/rnaseq_pipeline/tests/test_data/database_files/bioSample/bioSample_J.PLAGGENBERG_01.06.20.xlsx')
        known_col_dict = {'bioSampleNumber': 'biosampleNumber'}
        known_row_dict = {}
        self.assertDictEqual(known_col_dict, func_col_dict)
        self.assertDictEqual(known_row_dict, func_row_dict)

    def test_subdirectoryReport(self):
        db = DatabaseAccuracyObject(config_file='/home/chase/code/brentlab/rnaseq_pipeline/config/test_rnaseq_pipeline_config.ini')
        db.subdirectoryReport('bioSample', db, db.database_dict['bioSample'])
        self.assertTrue(os.path.isfile())

    def test_fullReport(self):
        dba = DatabaseAccuracyObject(config_file='/home/chase/code/brentlab/rnaseq_pipeline/config/test_rnaseq_pipeline_config.ini')
        dba.fullReport()

    def test_metadataAccuracyObject(self):
        raise NotImplementedError


if __name__ == '__main__':
    unittest.main()
