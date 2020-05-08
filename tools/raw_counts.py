#!/usr/bin/env python
import sys
import argparse
import os.path
import numpy as np
import glob
from rnaseq_tools import utils

# TODO: Update with OrganismData object (no more need to input gene list). Better commeting and explanation of each step
#  in createCountMatrix. Try to abstract this for use when making any matrix (counts qual_assess_1 and 2, etc)

def main(argv):
    args = parseArgs(argv)
#    if not os.path.exists(args.experiment_directory):
#        sys.exit('ERROR: %s does not exist.' % args.experiment_directory)
#    if not os.path.exists(args.gene_list):
#        sys.exit('ERROR: %s does not exist.' % args.gene_list)
#
#    count_matrix = createCountMatrix(args.experiment_directory, args.gene_list)
#
    exp_name = utils.dirName(args.experiment_directory)

    count_file_path = os.path.join(args.experiment_directory, exp_name + '_raw_count.csv')
    print(utils.dirName(args.experiment_directory))
    print(args.experiment_directory)
    print(count_file_path)
#
#    np.savetxt(count_file_path, count_matrix, delimiter=',', fmt='%s')

#    print(count_matrix[1:10])

def parseArgs(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--experiment_directory', required=True,
                        help='The directory created by create_experiment. The raw_count csv will be output in this directory with the name <experiment_dir/experiment_dirname>_raw_count.csv')
    parser.add_argument('-l', '--gene_list', required=True,
                        help='Gene list. This gene list should be a subset or the same set of annotated genes in GTF/GFF.')
    return parser.parse_args(argv[1:])

def createCountMatrix(exp_dir, gene_list):
    """
        create count matrix with list of genes as rows and samples as columns
        :param exp_dir: experiment directory created by create_experiment.py
        :param gene_list: see the OrganismData_config.ini in genome_files/<organism>
        :returns: A count matrix of all genes (rows) by all samples (columns)
    """
    exp_dir = utils.addForwardSlash(exp_dir)
    search_pattern = exp_dir + '*_read_count.tsv'
    sample_counts_list = glob.glob(search_pattern)
    # TODO: sanity check against num rows in query?

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
        # TODO: ask yiming about this line
        indx = np.where([gene_list[i] == count[:, 0] for i in range(len(gene_list))])[1]
        count_mtx = np.hstack((count_mtx, count[indx, 1].reshape(-1, 1)))
        sample_header = os.path.basename(file)
        header.append(sample_header)

    # combine count matrix and column headings
    count_mtx = np.vstack((header, count_mtx))

    return count_mtx


if __name__ == '__main__':
    main(sys.argv)
