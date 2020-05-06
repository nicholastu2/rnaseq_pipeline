import unittest
from unittest.mock import patch
from rnaseq_tools.OrganismData import OrganismData

class MyTestCase(unittest.TestCase):

    def setUp(self):
        self.patcher = patch('os.path.exists', return_value = True)
        self.patcher.start()

    @patch('rnaseq_tools.OrganismData')
    def test_standardDirectoryStructure(self, mock_organism_data):
        mock_organism_data.lts_rnaseq_data = '/lts/mblab/Crypto/rnaseq_data'
        mock_organism_data.pipeline_version = 'v1.0'
        mock_organism_data.mblab_scratch = '/scratch/mblab'
        mock_organism_data.scratch_database_files = '/scratch/mblab/mblab.shared/database_files'
        mock_organism_data.mblab_shared = '/scratch/mblab/mblab.shared'
        mock_organism_data._user = 'username'
        mock_organism_data.year_month_day = '20200506'
        mock_organism_data.logger = None
        mock_organism_data.standardDataStructure = OrganismData.standardDirectoryStructure
        mock_organism_data.patch('utils.executeSubProcess')
        mock_organism_data.patch('utils.mkdirp', return_value=True)
        #mock_organism_data.standardDataStructure(mock_organism_data)

        attributes = ['lts_rnaseq_data', 'pipeline_version', 'mblab_scratch', 'scratch_database_files',
                      'mblab_shared', 'lts_sequence', 'lts_align_expr', 'scratch_sequence', 'user_rnaseq_pipeline_directory',
                      'genome_files', 'reports', 'sbatch_log', 'log_dir', 'log_file', 'job_scripts', 'rnaseq_tmp',
                      'config_file', 'align_count_path']

        # check to see if all StandardData attributes have been set
        for assigned_attributes in attributes:
            self.assertTrue(hasattr(mock_organism_data, assigned_attributes))

        # check for OrganismData specific attributes
        org_specific = ['genome_files', 'organism_config_file', 'organism', 'feature_type']
        for org_attributes in org_specific:
            self.assertTrue(hasattr(mock_organism_data, org_attributes))


def tearDown(self):
    self.patcher.stop()


if __name__ == '__main__':
    unittest.main()