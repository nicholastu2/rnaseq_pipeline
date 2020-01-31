#!/usr/bin/env python
import sys
import argparse
import os.path
import numpy as np
import re
import glob
from utils import addForwardSlash


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--experiment_directory', required=True,
                        help='The directory created by create_experiment')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory for count matrix.')
    parser.add_argument('-l', '--gene_list', required=True,
                        help='Gene list. This gene list should be a subset or the same set of annotated genes in GTF/GFF.')
    return parser.parse_args(argv[1:])

def createCountMatrix(exp_dir, gene_list):
    # create count matrix with list of genes as rows and samples as columns
    # Args: sample_summary.csv (see queryDB in rnaseq_pipeline/tools)
    # Returns: A count matrix of all genes (rows) by all samples (columns)

    exp_dir = addForwardSlash(exp_dir)
    search_pattern = exp_dir + '*_read_count.tsv'
    sample_counts_list = glob.glob(search_pattern)
    # verify that the files exist
    for file_path in sample_counts_list:
        if not os.path.exists(file_path):
            sys.exit('ERROR: %s does not exist.\n... Aborted preparing count matrix.' % file_path)

    # import list of genes from genome
    gene_list = np.loadtxt(gene_list, dtype=str)

    # create count_mtx with rownames gene_list
    count_mtx = gene_list.reshape(-1, 1)
    # create list of column header names
    header = ['gene_id']

    # append samples as columns to count_mtx 
    for file in sample_counts_list:
        # append sample column counts
        print('... working on %s' % file)
        count = np.loadtxt(file, dtype=str)
        indx = np.where([gene_list[i] == count[:, 0] for i in range(len(gene_list))])[1]
        count_mtx = np.hstack((count_mtx, count[indx, 1].reshape(-1, 1)))
        # append sample column heading name 
        sample_header = re.sub(".fastq.gz|.fastq", "", file)
        header.append(sample_header)

    # combine count matrix and column headings
    count_mtx = np.vstack((header, count_mtx))

    return count_mtx

def main(argv):
    parsed = parse_args(argv)
    if not os.path.exists(parsed.experiment_directory):
        sys.exit('ERROR: %s does not exist.' % parsed.input_lookup)
    if not os.path.exists(parsed.gene_list):
        sys.exit('ERROR: %s does not exist.' % parsed.gene_list)

    count_matrix = createCountMatrix(parsed.experiment_directory, parsed.gene_list)

    count_file_path = os.path.join(parsed.output, os.path.basename(parsed.experiment_directory)) + 'gene_count.csv'

    np.savetxt(count_file_path, count_matrix, delimiter=',', fmt='%s')

    print(count_matrix[1:10])


if __name__ == '__main__':
    main(sys.argv)
