import utils
import pandas as pd
import subprocess

def makeIgvSnapshotDict(sample_list, query_df, wildtype):
    igv_snapshot_dict = {}
    genotype_list = []
    bam_file_list = []
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

query_df = pd.read_csv('/home/chase/Documents/CNAG_05420_all_after_update.csv')

wildtype = 'CNAG_05420'

cmd = 'ml rnaseq_pipeline'
utils.executeSubProcess(cmd)

print(makeIgvSnapshotDict(sample_list, query_df, wildtype))

# >>> class TestConfigure:
# ...     def __init__(self, config_file):
# ...         self.config_file = config_file
# ...         self.self_type = 'StandardData'
# ...         utils.configure(self)
