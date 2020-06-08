import unittest
from rnaseq_tools import annotation_tools

class MyTestCase(unittest.TestCase):
    def test_parseGtf_new(self):

        known_dict = {'CNAG_00001': {'chromosome': 'chr1', 'strand': '-',
                                     'biotype': 'hypothetical protein',
                                     'gc_content': 179, 'gene_start': 11499,
                                     'gene_stop': 11065, 'cds_length': 327},

                      'CNAG_00002': {'chromosome': 'chr1', 'strand': '+', 'biotype': 'hypothetical protein', 'gc_content': 249, 'gene_start': 12800, 'gene_stop': 13581, 'cds_length': 525}}

        parsed_dict = annotation_tools.parseGtf_new('/home/chase/code/brentlab/rnaseq_pipeline/tests/test_data/KN99_2_genes.gtf', '/home/chase/Desktop/rnaseq_pipeline/rnaseq_pipeline/genome_files/KN99/crNeoKN99.fasta')

        self.assertDictEqual(known_dict, parsed_dict)


if __name__ == '__main__':
    unittest.main()
