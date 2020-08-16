import unittest
from unittest.mock import patch
from rnaseq_tools.OrganismDataObject import OrganismData
import glob

class MyTestCase(unittest.TestCase):
    def test_createCountSheet(self):
        config_file = '/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline_config.ini'

        od = OrganismData(config_file=config_file, interactive=True, organism='S288C_R64')

        path='/home/chase/Desktop/rnaseq_pipeline/run_4027_samples/count/*_read_count.tsv'
        count_file_list = glob.glob(path)

        df = od.createCountSheet(count_file_list=count_file_list)
        print(df)


if __name__ == '__main__':
    unittest.main()