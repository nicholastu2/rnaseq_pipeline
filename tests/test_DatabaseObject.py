import unittest
import os
import pandas as pd
import logging
import sys
from unittest.mock import patch
from rnaseq_tools import utils
from rnaseq_tools.DatabaseObject import DatabaseObject


class MyTestCase(unittest.TestCase):
    def test_DataObjectConstructor(self):
        db_object = DatabaseObject('/path/to/database/top')
        self.assertEqual(db_object.database_top, '/path/to/database/top')

    def test_setDataDirDict(self):
        known_file = ['database-files/bioSample/bioSample_J.PLAGGENBERG_07.24.19.xlsx', 'database-files/bioSample/bioSample_hbrown_09.10.19.xlsx']
        db_object = DatabaseObject('/path/to/database/top', database_subdirectories=['known_key'])
        db_object.setDatabaseDict()
        self.assertDictEqual({'known_key': known_file}, db_object.database_dict)

    def test_setConcatDatabaseDict(self):
        db_object = DatabaseObject(config_file='/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline_config.ini', interactive=True)
        db_object.setDatabaseDict()
        db_object.setConcatDatabaseDict()
        known = {'col_1': [3,2,3,2], 'col_2': ['a','b','a','b']}
        df_dict = db_object.concat_database_dict['known_key'].to_dict('list')
        self.assertDictEqual(known, df_dict)


    def test_setDatabaseDataframe(self):
        db_object = DatabaseObject('/path/to/database/top', database_subdirectories=['known_1', 'known_2', 'known_3'])
        db_object.setDatabaseDict()
        db_object.setDatabaseDataframe()
        col1_lst = [[[3]*4, [2]*4]*2]
        col1_lst = [item for sublist in col1_lst for item in sublist]
        col1_lst = [item for sublist in col1_lst for item in sublist]
        col2_lst = [[['a']*4, ['b']*4]*2]
        col2_lst = [item for sublist in col2_lst for item in sublist]
        col2_lst = [item for sublist in col2_lst for item in sublist]
        known_dict = {'col_1': col1_lst,'col_2': col2_lst}
        df_dict = db_object.database_df.to_dict('list')
        self.assertDictEqual(known_dict, df_dict)

    def test_setKeyColumns(self):
        db_object = DatabaseObject('/path/to/database/top', database_subdirectories=['known_1', 'known_2', 'known_3'])
        db_object.setDatabaseDict()
        db_object.setDatabaseDataframe()
        known_keys = [['col_1', 'col_2'],['col_1', 'col_2']]
        db_object_key_columns = [list(db_object.database_key_columns[0]), list(db_object.database_key_columns[1])]
        self.assertListEqual(known_keys, db_object_key_columns)

    def test_setFilterJson(self):
        pass

    def test_dropRowsIfEmptyFastqFilename(self):
        pass

    def test_filterDatabaseDataframe(self):
        pass


if __name__ == '__main__':
    unittest.main()
