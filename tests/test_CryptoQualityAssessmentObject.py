import unittest
from rnaseq_tools.CryptoQualityAssessmentObject import CryptoQualityAssessmentObject
from rnaseq_tools import utils
from glob import glob

class MyTestCase(unittest.TestCase):
    def test_auditQualAssessDataFrame(self):
        qa = CryptoQualityAssessmentObject(config_file='/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline_config.ini',
                                           interactive=True)
        bam_file_list = utils.extractFiles('/home/chase/Desktop/tmp/run_1045_samples/align', '.bam')
        qual_asses_df_path = '/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline/experiments/run_1045_quality_summary.csv'

        query_df_path = '/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline/query/crypto_qa_tester_nextflow_full_query.csv'

        df = qa.auditQualAssessDataFrame(query_df_path, qual_asses_df_path, bam_file_list)
        print('yes')


if __name__ == '__main__':
    unittest.main()
