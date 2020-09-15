import unittest
from rnaseq_tools.IgvObject import IgvObject

class MyTestCase(unittest.TestCase):

    def test_igv_constructor(self):
        sample_sheet_path = '/home/chase/Desktop/tmp/igv_perturbed_query_tester.csv'
        wildtype_sheet_path = '/home/chase/Desktop/tmp/igv_wildtype_query_tester.csv'
        organism = 'KN99'
        scratch_alignment_source = '/mnt/htcf_scratch/chasem/rnaseq_pipeline/align_count_results'
        igv_output_dir = '/home/chase/Desktop/tmp/igv_output_test'
        igv_object = IgvObject(sample_sheet_path=sample_sheet_path,
                               wildtype_sheet_path=wildtype_sheet_path,
                               organism=organism,
                               scratch_alignment_source=scratch_alignment_source,
                               igv_output_dir=igv_output_dir,
                               config_file='/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline_config.ini',
                               interactive=True)
        wildtype_dict_from_sheet = {'37C.CO2_90': '1874_Brent_1_TGAGGTT_S1_R1_001.fastq.gz'}
        self.assertDictEqual(wildtype_dict_from_sheet, igv_object.wildtype_dict)

    def test_igv_constructor(self):
        sample_sheet_path = '/home/chase/Desktop/tmp/igv_perturbed_query_tester.csv'
        wildtype_sheet_path = '/home/chase/Desktop/tmp/igv_wildtype_query_tester.csv'
        organism = 'KN99'
        scratch_alignment_source = '/mnt/htcf_scratch/chasem/rnaseq_pipeline/align_count_results'
        igv_output_dir = '/home/chase/Desktop/tmp/igv_output_test'
        igv_object = IgvObject(sample_sheet_path=sample_sheet_path,
                               wildtype_sheet_path=wildtype_sheet_path,
                               organism=organism,
                               scratch_alignment_source=scratch_alignment_source,
                               igv_output_dir=igv_output_dir,
                               config_file='/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline_config.ini',
                               interactive=True)
        igv_object.makeIgvSnapshotDict()
        igv_object.writeIgvJobScript()
        print('yes')


if __name__ == '__main__':
    unittest.main()
