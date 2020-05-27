import unittest
from rnaseq_tools.QualityAssessmentObject import QualityAssessmentObject

class MyTestCase(unittest.TestCase):
    def test_parseAlignmentFile(self):
        x = QualityAssessmentObject.parseAlignmentLog('/mnt/htcf_scratch/chasem/rnaseq_pipeline/re_align_count/run_1028_KN99_no/run_1028_s_4_withindex_sequence_GGTCCTC_novoalign.log')
        print(x)

        self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()
