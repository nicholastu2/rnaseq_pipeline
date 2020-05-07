import unittest
import os
from unittest.mock import patch
from rnaseq_tools import utils
from rnaseq_tools.DatabaseObject import DatabaseObject

class MyTestCase(unittest.TestCase):
    def setUp(self):
        # for utils.extractTopmostFiles
        files = ['~database-files/bioSample/bioSample_hbrown_07.24.19.xlsx',
                 '._database-files/bioSample/bioSample_hbrown_09.13.19.xlsx',
                 '~&database-files/bioSample/bioSample_J.PLAGGENBERG_01.10.20.xlsx',
                 'database-files/bioSample/bioSample_J.PLAGGENBERG_07.24.19.xlsx']

        self.patcher = patch('os.path.exists', return_value=True)
        self.patcher = patch('rnaseq_tools.utils.extractTopmostFiles', return_value=files)
        self.patcher.start()

    def test_DataObjectConstructor(self):
        known_response = ['database-files/bioSample/bioSample_J.PLAGGENBERG_07.24.19.xlsx']
        db_object = DatabaseObject('/path/to/database/top')
        logger_path = db_object.logger_path
        self.assertTrue(os.path.isfile(logger_path))
        cmd = 'rm %s' %logger_path
        utils.executeSubProcess(cmd)
        self.assertEqual(db_object.database_top, '/path/to/database/top')

    def test_setDataDirDict(self):
        known_file = ['database-files/bioSample/bioSample_J.PLAGGENBERG_07.24.19.xlsx']
        db_object = DatabaseObject('/path/to/database/top', database_subdirectories = ['known_key'])
        db_object.setDatabaseDict()
        self.assertDictEqual({'known_key': known_file}, db_object.database_dict)

def tearDown(self):
    self.patcher.stop()


if __name__ == '__main__':
    unittest.main()
