import unittest
from rnaseq_tools.CryptoQualityAssessmentObject import CryptoQualityAssessmentObject
from rnaseq_tools.S288C_R64QualityAssessmentObject import S288C_R64QualityAssessmentObject
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
        crypto_qa = CryptoQualityAssessmentObject(bam_file_list=bam_list,
                                                  count_file_list=count_list,
                                                  novoalign_log_list=novoalign_logs,
                                                  coverage_check_flag=True,
                                                  query_path=query_path,
                                                  config_file='/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline_config.ini',
                                                  interactive=True)

        print('...compiling alignment information')
        # create dataframes storing the relevant alignment and count metadata from the novoalign and htseq logs
        crypto_qual_assess_1_df = crypto_qa.compileAlignCountMetadata()

        print('writing output to %s' % output_full_path)
        crypto_qual_assess_1_df.to_csv(output_full_path, index=False)

    def test_yeastQualAssess(self):

        output_path=''

        path_to_yeast_data_dir = '/home/chase/Desktop/rnaseq_pipeline/crypto_test_data/*/*/*.'

        bam_list = glob.glob(path_to_yeast_data_dir+'bam', recursive=True)
        count_list = glob.glob(path_to_yeast_data_dir+'tsv', recursive=True)
        novoalign_logs = glob.glob(path_to_yeast_data_dir[:-1]+'novoalign.log', recursive=True)
        query_path = '/home/chase/code/brentlab/rnaseq_pipeline/tests/test_data/yeast_qual_assess_test_query.csv'

        crypto_qa = S288C_R64QualityAssessmentObject(bam_file_list=bam_list,
                                                     count_file_list=count_list,
                                                     novoalign_log_list=novoalign_logs,
                                                     coverage_check_flag=True,
                                                     query_path=query_path,
                                                     config_file='/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline_config.ini',
                                                     interactive=True)

        print('...compiling alignment information')
        # create dataframes storing the relevant alignment and count metadata from the novoalign and htseq logs
        crypto_qual_assess_1_df = crypto_qa.compileAlignCountMetadata()

        print('writing output to %s' % output_path)
        crypto_qual_assess_1_df.to_csv(output_path, index=False)


    def test_auditQualAssessDataFrame(self):
        qa = CryptoQualityAssessmentObject(config_file='/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline_config.ini',
                                           interactive=True)
        bam_file_list = utils.extractFiles('/home/chase/Desktop/tmp/run_1045_samples/align', '.bam')
        qual_asses_df_path = '/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline/experiments/run_1045_quality_summary.csv'

        query_df_path = '/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline/query/combined_df_20200720.csv'

        df = qa.auditQualAssessDataFrame(query_df_path, qual_asses_df_path, bam_file_list)
        print('yes')


if __name__ == '__main__':
    unittest.main()
