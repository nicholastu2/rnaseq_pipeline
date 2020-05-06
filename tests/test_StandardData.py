import unittest
from unittest.mock import patch
from rnaseq_tools.StandardData import StandardData

class MyTestCase(unittest.TestCase):

    def setUp(self):
        self.patcher = patch('os.path.exists', return_value = True)
        self.patcher.start()

    @patch('rnaseq_tools.StandardData')
    def test_standardDirectoryStructure(self, mock_standard_data):
        mock_standard_data.lts_rnaseq_data = '/lts/mblab/Crypto/rnaseq_data'
        mock_standard_data.pipeline_version = 'v1.0'
        mock_standard_data.mblab_scratch = '/scratch/mblab'
        mock_standard_data.scratch_database_files = '/scratch/mblab/mblab.shared/database_files'
        mock_standard_data.mblab_shared = '/scratch/mblab/mblab.shared'
        mock_standard_data._user = 'username'
        mock_standard_data.year_month_day = '20200506'
        mock_standard_data.logger = None
        mock_standard_data.standardDataStructure = StandardData.standardDirectoryStructure
        mock_standard_data.patch('utils.executeSubProcess')
        mock_standard_data.patch('utils.mkdirp', return_value=True)
        mock_standard_data.standardDataStructure(mock_standard_data)

        attributes = ['lts_rnaseq_data', 'pipeline_version', 'mblab_scratch', 'scratch_database_files',
                      'mblab_shared', 'lts_sequence', 'lts_align_expr', 'scratch_sequence', 'user_rnaseq_pipeline_directory',
                      'genome_files', 'reports', 'sbatch_log', 'log_dir', 'log_file', 'job_scripts', 'rnaseq_tmp',
                      'config_file', 'align_count_path']

        # check to see if all attributes have been set
        for assigned_attributes in attributes:
            self.assertTrue(hasattr(mock_standard_data, assigned_attributes))


def tearDown(self):
    self.patcher.stop()


if __name__ == '__main__':
    unittest.main()
