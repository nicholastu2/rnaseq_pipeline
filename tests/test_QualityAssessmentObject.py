import unittest
from rnaseq_tools.QualityAssessmentObject import QualityAssessmentObject


class MyTestCase(unittest.TestCase):
    def test_coverageCheck(self):
        nextflow_list_of_files = ['blah_htseq_log', 'sequence/run_0673_samples/run_673_s_4_withindex_sequence_TGAGGTT_sorted_aligned_reads.bam',
                                  'sequence/run_0673_samples/run_673_s_4_withindex_sequence_GCTTAGA_sorted_aligned_reads.bam',
                                  'Brent_large_3b-1_GTAC_1_TGAGGTT_universal_S38_R1_001_sorted_aligned_reads.bam',
                                  'Brent_large_3b-2_GTAC_2_GCTTAGA_universal_S39_R1_001_sorted_aligned_reads.bam',
                                  'Brent_01_GTAC1_SIC_Index10_TGAGGTTATC_GCTTCTGT_S2_R1_001_sorted_aligned_reads.bam',
                                  'Brent_02_GTAC2_SIC_Index10_GCTTAGAATC_GCTTCTGT_S41_R1_001_sorted_aligned_reads.bam',
                                  'sequence/run_0759_samples/run_759_s_1_withindex_sequence_AGTTATG_sorted_aligned_reads.bam',
                                  'sequence/run_0648_samples/s_7_withindex_sequence_GCCTCCC_sorted_aligned_reads.bam']

        qa = QualityAssessmentObject(query_path='/home/chase/code/brentlab/rnaseq_pipeline/tests/test_user/rnaseq_pipeline/query/test_database_df.csv',
                                     coverage_check_flag=True,
                                     nextflow_list_of_files=nextflow_list_of_files,
                                     config_file='/home/chase/code/brentlab/rnaseq_pipeline/config/test_rnaseq_pipeline_config.ini',
                                     interactive=True)
        bam_file_paths = [bam_file for bam_file in qa.nextflow_list_of_files if '_sorted_aligned_reads.bam' in bam_file]
        qa.coverageCheck()
        print(bam_file_paths)

    def test_igvShot(self):
        qa = QualityAssessmentObject(config_file='/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline_config.ini',
                                     interactive=True)
        bam_path_1='/mnt/htcf_scratch/chasem/rnaseq_pipeline/align_count_results/run_1658_samples/align/1658_Brent_Crypto_TimeCourse2016_26_TCAACTG_S25_R1_001_sorted_aligned_reads_with_annote.bam'
        bam_path_2='/mnt/htcf_scratch/chasem/rnaseq_pipeline/align_count_results/run_1658_samples/align/1658_Brent_Crypto_TimeCourse2016_27_TGTTTGT_S26_R1_001_sorted_aligned_reads_with_annote.bam'
        bam_list=[bam_path_1, bam_path_2]
        gene_list=['CNAG_NAT', 'CNAG_G418']
        qa.igvShot(bam_list, organism='KN99', gene=gene_list)


if __name__ == '__main__':
    unittest.main()
