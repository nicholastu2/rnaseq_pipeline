from rnaseq_tools import utils
import pandas as pd
import subprocess

def makeIgvSnapshotDict(sample_list, query_df, wildtype):
    igv_snapshot_dict = {}
    genotype_list = []
    wildtype_sample_list = []
    for sample in sample_list:
        # if sample.endswith('.fastq.gz'):
        #     sample = utils.pathBaseName(sample) + '_read_count.tsv'
        query_row_with_sample = query_df[query_df['fastqFileName'] == sample]
        # extract value in genotype column
        genotype = query_row_with_sample['genotype'].values[0]
        # split on period if this is a double perturbation. Regardless of whether a . is present, genotype will be cast to a list eg ['CNAG_00000'] or ['CNAG_05420', 'CNAG_01438']
        genotype = genotype.split('.')
        # store the wt sample, and then move onto the next sample in the list
        if not genotype[0] == wildtype:
            genotype_list.extend(genotype)
            # create bamfile name
            bamfile = sample.replace('_read_count.tsv', '_sorted_aligned_reads.bam')
            # add to igv_snapshot_dict
            igv_snapshot_dict.setdefault(sample, {}).setdefault('gene', []).extend(genotype)
            igv_snapshot_dict[sample]['bam'] = bamfile
            igv_snapshot_dict[sample]['bed'] = None
        else:
            wt_sample = sample
    # if the wt genotype was found, create entry in the following form {sample_read_counts.tsv: {'gene': [perturbed_gene_1, perturbed_gene_2, ...], 'bam': wt.bam, 'bed': created_bed.bed}
    if wt_sample:
        igv_snapshot_dict.setdefault(wt_sample, {}).setdefault('gene', []).extend(genotype_list)
        igv_snapshot_dict[wt_sample]['bam'] = bamfile
        igv_snapshot_dict[wt_sample]['bed'] = None
    return igv_snapshot_dict

sample_list = ['sequence/run_1045_samples/run_1045_s_7_withindex_sequence_CACCTCC.fastq.gz',
'sequence/run_1045_samples/run_1045_s_7_withindex_sequence_ATCGAGC.fastq.gz',
'sequence/run_1045_samples/run_1045_s_7_withindex_sequence_TACTCTA.fastq.gz',
'sequence/run_1045_samples/run_1045_s_7_withindex_sequence_AGACTGA.fastq.gz',
'sequence/run_1045_samples/run_1045_s_7_withindex_sequence_CTTGGAA.fastq.gz',
]

sample_list = [
'sequence/run_1800_samples/1800_Brent_12_TAACAAG_S12_R1_001.fastq.gz',
'sequence/run_1800_samples/1800_Brent_13_GAGGCGT_S13_R1_001.fastq.gz',
'sequence/run_1800_samples/1800_Brent_14_TTTAACT_S14_R1_001.fastq.gz',
'sequence/run_1800_samples/1800_Brent_15_GGTCCTC_S15_R1_001.fastq.gz',
'sequence/run_1800_samples/1800_Brent_19_GAGTACG_S19_R1_001.fastq.gz',
'sequence/run_1800_samples/1800_Brent_20_ACAGATA_S20_R1_001.fastq.gz',
'sequence/run_1800_samples/1800_Brent_21_CTCAATG_S21_R1_001.fastq.gz',
'sequence/run_0773_samples/run_773_s_2_withindex_sequence_GGCAGCG.fastq.gz',
'sequence/run_0773_samples/run_773_s_2_withindex_sequence_CCATCAT.fastq.gz',
'sequence/run_0773_samples/run_773_s_2_withindex_sequence_TAACAAG.fastq.gz',
'sequence/run_0773_samples/run_773_s_2_withindex_sequence_GAGTACG.fastq.gz',
]

from rnaseq_tools.IgvObject import IgvObject
sample_list = [
'sequence/run_1800_samples/1800_Brent_12_TAACAAG_S12_R1_001.fastq.gz',
'sequence/run_1800_samples/1800_Brent_13_GAGGCGT_S13_R1_001.fastq.gz',
'sequence/run_1800_samples/1800_Brent_14_TTTAACT_S14_R1_001.fastq.gz',
'sequence/run_1800_samples/1800_Brent_15_GGTCCTC_S15_R1_001.fastq.gz',
'sequence/run_1800_samples/1800_Brent_19_GAGTACG_S19_R1_001.fastq.gz',
'sequence/run_1800_samples/1800_Brent_20_ACAGATA_S20_R1_001.fastq.gz',
'sequence/run_1800_samples/1800_Brent_21_CTCAATG_S21_R1_001.fastq.gz',
'sequence/run_0773_samples/run_773_s_2_withindex_sequence_GGCAGCG.fastq.gz',
'sequence/run_0773_samples/run_773_s_2_withindex_sequence_CCATCAT.fastq.gz',
'sequence/run_0773_samples/run_773_s_2_withindex_sequence_TAACAAG.fastq.gz',
'sequence/run_0773_samples/run_773_s_2_withindex_sequence_GAGTACG.fastq.gz',
]

query_path = '/scratch/mblab/chasem/old_rnaseq/query/CNAG_05420_all.csv'

wildtype = 'CNAG_05420'

# expirement directory already make -- see constructor
igv = IgvObject(query_sheet_path = query_path, sample_list = sample_list,wildtype = 'CNAG_05420',experiment_dir = '/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/20200322_120638_igv_files', organism = 'KN99', email = 'chasem@wustl.edu')

query_df = pd.read_csv('/home/chase/Documents/CNAG_05420_all_after_update.csv')


print(makeIgvSnapshotDict(sample_list, query_df, wildtype))


# sdf = StandardData(query_sheet_path = query_df)
#
# def extractValueFromStandardRow(self, filter_column, filter_value, extract_column, run_num_with_leading_zero=False):
# def extractValueFromStandardRow(self, filter_column, filter_value, extract_column, run_num_with_leading_zero=False):
#     """
#     extract a value from a row (selected by filter_value) of self.query_df
#     :param filter_column:
#     :param filter_value:
#     :param extract_column:
#     :param run_num_with_leading_zero: if true, return the run number as a string with a leading 0 if it is in self._run_numbers_with_zeros
#     :returns: a value extracted from a certain column of a certain row
#     """
#     row = self.query_df[self.query_df[filter_column] == filter_value]
#
#     extracted_value = row[extract_column].values(0)
#
#     if run_num_with_leading_zero:
#         if extracted_value in self._run_numbers_with_zeros:
#             extracted_value = self._run_numbers_with_zeros[extracted_value]
#
#     return extracted_value


{'run_1045_s_7_withindex_sequence_CACCTCC_read_count.tsv':
       {'gene': ['CNAG_05420_over', 'CNAG_05420_over'],
        'bam': '/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/20200320_161653_igv_files/run_1045_s_7_withindex_sequence_CACCTCC_sorted_aligned_reads.bam',
        'bed': None},
 'run_1045_s_7_withindex_sequence_ATCGAGC_read_count.tsv':
       {'gene': ['CNAG_05420_over', 'CNAG_05420_over'],
        'bam': '/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/20200320_161653_igv_files/run_1045_s_7_withindex_sequence_ATCGAGC_sorted_aligned_reads.bam',
        'bed': None},
'run_1045_s_7_withindex_sequence_TACTCTA_read_count.tsv':
       {'gene': ['CNAG_05420_over', 'CNAG_05420_over'],
        'bam': '/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/20200320_161653_igv_files/run_1045_s_7_withindex_sequence_TACTCTA_sorted_aligned_reads.bam',
        'bed': None},
'run_1045_s_7_withindex_sequence_AGACTGA_read_count.tsv':
       {'gene': ['CNAG_05420_over'],
        'bam': '/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/20200320_161653_igv_files/run_1045_s_7_withindex_sequence_AGACTGA_sorted_aligned_reads.bam',
        'bed': None},
'run_1045_s_7_withindex_sequence_CTTGGAA_read_count.tsv': {'gene': ['CNAG_05420_over'],
        'bam': '/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/20200320_161653_igv_files/run_1045_s_7_withindex_sequence_CTTGGAA_sorted_aligned_reads.bam',
        'bed': None}}

{'1800_Brent_12_TAACAAG_S12_R1_001_read_count.tsv':
     {'gene': ['CNAG_04864', 'CNAG_05420', 'CNAG_04864', 'CNAG_05420', 'CNAG_04864', 'CNAG_05420', 'CNAG_01438', 'CNAG_05420', 'CNAG_01438', 'CNAG_05420', 'CNAG_01438', 'CNAG_05420', 'CNAG_01551', 'CNAG_05420', 'CNAG_01551', 'CNAG_05420', 'CNAG_01551', 'CNAG_05420', 'CNAG_05222', 'CNAG_05420'],
      'bam': '/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/20200320_173101_igv_files/1800_Brent_12_TAACAAG_S12_R1_001_sorted_aligned_reads.bam',
      'bed': None},
 '1800_Brent_13_GAGGCGT_S13_R1_001_read_count.tsv':
     {'gene': ['CNAG_04864', 'CNAG_05420'],
      'bam': '/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/20200320_173101_igv_files/1800_Brent_13_GAGGCGT_S13_R1_001_sorted_aligned_reads.bam',
      'bed': None},
 '1800_Brent_14_TTTAACT_S14_R1_001_read_count.tsv':
     {'gene': ['CNAG_04864', 'CNAG_05420'],
      'bam': '/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/20200320_173101_igv_files/1800_Brent_14_TTTAACT_S14_R1_001_sorted_aligned_reads.bam',
      'bed': None}, '1800_Brent_15_GGTCCTC_S15_R1_001_read_count.tsv': {'gene': ['CNAG_04864', 'CNAG_05420'], 'bam': '/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/20200320_173101_igv_files/1800_Brent_15_GGTCCTC_S15_R1_001_sorted_aligned_reads.bam', 'bed': None}, '1800_Brent_19_GAGTACG_S19_R1_001_read_count.tsv': {'gene': ['CNAG_01438', 'CNAG_05420'], 'bam': '/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/20200320_173101_igv_files/1800_Brent_19_GAGTACG_S19_R1_001_sorted_aligned_reads.bam', 'bed': None}, '1800_Brent_20_ACAGATA_S20_R1_001_read_count.tsv': {'gene': ['CNAG_01438', 'CNAG_05420'], 'bam': '/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/20200320_173101_igv_files/1800_Brent_20_ACAGATA_S20_R1_001_sorted_aligned_reads.bam', 'bed': None}, '1800_Brent_21_CTCAATG_S21_R1_001_read_count.tsv': {'gene': ['CNAG_01438', 'CNAG_05420'], 'bam': '/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/20200320_173101_igv_files/1800_Brent_21_CTCAATG_S21_R1_001_sorted_aligned_reads.bam', 'bed': None}, 'run_773_s_2_withindex_sequence_GGCAGCG_read_count.tsv': {'gene': ['CNAG_01551', 'CNAG_05420'], 'bam': '/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/20200320_173101_igv_files/run_773_s_2_withindex_sequence_GGCAGCG_sorted_aligned_reads.bam', 'bed': None}, 'run_773_s_2_withindex_sequence_CCATCAT_read_count.tsv': {'gene': ['CNAG_01551', 'CNAG_05420'], 'bam': '/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/20200320_173101_igv_files/run_773_s_2_withindex_sequence_CCATCAT_sorted_aligned_reads.bam', 'bed': None}, 'run_773_s_2_withindex_sequence_TAACAAG_read_count.tsv': {'gene': ['CNAG_01551', 'CNAG_05420'], 'bam': '/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/20200320_173101_igv_files/run_773_s_2_withindex_sequence_TAACAAG_sorted_aligned_reads.bam', 'bed': None}, 'run_773_s_2_withindex_sequence_GAGTACG_read_count.tsv': {'gene': ['CNAG_05222', 'CNAG_05420'], 'bam': '/scratch/mblab/chasem/rnaseq_pipeline/rnaseq_tmp/20200320_173101_igv_files/run_773_s_2_withindex_sequence_GAGTACG_sorted_aligned_reads.bam', 'bed': None}}

