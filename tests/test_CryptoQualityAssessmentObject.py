import unittest
from rnaseq_tools.CryptoQualAssessAuditObject import CryptoQualAssessAuditObject
from rnaseq_tools.CryptoQualityAssessmentObject import CryptoQualityAssessmentObject
from rnaseq_tools.S288C_R54QualAssessAuditObject import S288C_R54QualAssessAuditObject
from rnaseq_tools import utils
from glob import glob

class MyTestCase(unittest.TestCase):
    def test_cryptoQualAssess(self):

        output_full_path='/home/chase/Desktop/rnaseq_pipeline/crypto_test_qa_1.csv'

        path_to_crypto_data_dir = '/home/chase/Desktop/rnaseq_pipeline/crypto_test_data/*/*/*.'

        bam_list = glob(path_to_crypto_data_dir+'bam', recursive=True)
        count_list = glob(path_to_crypto_data_dir+'tsv', recursive=True)
        novoalign_logs = glob(path_to_crypto_data_dir[:-1]+'novoalign.log', recursive=True)
        query_path = '/home/chase/code/brentlab/rnaseq_pipeline/tests/test_data/crypto_qual_assess_test_query.csv'
        crypto_qa = CryptoQualityAssessmentObject(organism='KN99',
                                                  bam_file_list=bam_list,
                                                  count_file_list=count_list,
                                                  novoalign_log_list=novoalign_logs,
                                                  coverage_check_flag=True,
                                                  query_path=query_path,
                                                  config_file='/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline_config.ini',
                                                  interactive=True)



        print('writing output to %s' % output_full_path)
        crypto_qa.qual_assess_df.to_csv(output_full_path, index=False)

    def test_CryptoQualAssessAuditObject(self):

        output_full_path='/home/chase/Desktop/rnaseq_pipeline/crypto_test_qa_1.csv'

        path_to_crypto_data_dir = '/home/chase/Desktop/rnaseq_pipeline/crypto_test_data/*/*/*.'

        bam_list = glob(path_to_crypto_data_dir+'bam', recursive=True)
        count_list = glob(path_to_crypto_data_dir+'tsv', recursive=True)
        novoalign_logs = glob(path_to_crypto_data_dir[:-1]+'novoalign.log', recursive=True)
        query_path = '/home/chase/code/brentlab/rnaseq_pipeline/tests/test_data/crypto_qual_assess_test_query.csv'
        crypto_audit = CryptoQualAssessAuditObject(organism='KN99',
                                                  bam_file_list=bam_list,
                                                  count_file_list=count_list,
                                                  novoalign_log_list=novoalign_logs,
                                                  coverage_check_flag=True,
                                                  query_path=query_path,
                                                  config_file='/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline_config.ini',
                                                  interactive=True)



        print('writing output to %s' % output_full_path)
        crypto_audit.qual_assess_df.to_csv(output_full_path, index=False)

    def test_yeastQualAssess(self):

        output_path='/home/chase/Desktop/rnaseq_pipeline/yeast_test_qa_1.csv'

        path_to_yeast_data_dir = '/home/chase/Desktop/rnaseq_pipeline/run_4027_samples/*/*.'

        bam_list = glob(path_to_yeast_data_dir+'bam', recursive=True)
        count_list = glob(path_to_yeast_data_dir+'tsv', recursive=True)
        novoalign_logs = glob(path_to_yeast_data_dir[:-1]+'novoalign.log', recursive=True)
        query_path = '/home/chase/Desktop/rnaseq_pipeline/yeast_tester.csv'

        yeast_qa = S288C_R54QualAssessAuditObject(organism='S288C_R64',
                                                     bam_file_list=bam_list,
                                                     count_file_list=count_list,
                                                     novoalign_log_list=novoalign_logs,
                                                     coverage_check_flag=False,
                                                     query_path=query_path,
                                                     config_file='/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline_config.ini',
                                                     interactive=True)

        print('...compiling alignment information')

        print('writing output to %s' % output_path)
        yeast_qa.qual_assess_df.to_csv(output_path, index=False)


if __name__ == '__main__':
    unittest.main()
