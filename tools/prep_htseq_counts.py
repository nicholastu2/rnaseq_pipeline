#!/usr/bin/env python
import sys
import argparse
import os.path
import pandas as pd
import numpy as np
import re

def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query', required=True,
					help='query the database using queryDB.py for the libraries you wish to analyze')
    parser.add_argument('-o', '--output_count_file', required=True,
					help='Output count matrix.')
    parser.add_argument('-l', '--gene_list', required=True,
						help='Gene list. This gene list should be a subset or the same set of annotated genes in GTF/GFF.')
    return parser.parse_args(argv[1:])


def getFastqBasename(file):

    # create list of samples from the fastqFileName column of sample_summary.csv
    # Args: sample_summary.csv (see queryDB in rnaseq_pipeline/tools)
    # Returns: a list of the samples (the fastq file names minus the file extentions)

    # extract column fastqFileName of sample_summary.csv as list
    sample_fastq_list = pd.read_csv(file)['fastqFileName']
    
    # eliminate file extension (either .fastq.gz or .fastq) 
    sample_basename = [re.sub(".fastq.gz|.fastq", "", basename) for basename in sample_fastq_list]
    
    return sample_basename

def create_count_matrix(query, gene_list):
	
    # create count matrix with list of genes as rows and samples as columns
    # Args: sample_summary.csv (see queryDB in rnaseq_pipeline/tools)
    # Returns: A count matrix of all genes (rows) by all samples (columns)
    
    # get list of sample file basenames from fastq filenames
    sample_basename = getFastqBasename(query)
    # concatenate on filepath and extention
    sample_counts_list = ["alignment/novoalign/"+ basename +"_read_count.tsv" for basename in sample_basename]
    # verify that the files exist
    for file_path in sample_counts_list:
        if not os.path.exists(file_path):
            sys.exit('ERROR: %s does not exist.\n... Aborted preparing count matrix.' % file_path)
        
    # import list of genes from genome
    gene_list = np.loadtxt(gene_list, dtype=str)

    # create count_mtx with rownames gene_list
    count_mtx = gene_list.reshape(-1,1)
    # create list of column header names
    header = ['gene_id']

    # append samples as columns to count_mtx 
    for file in sample_counts_list:
        # append sample column counts
        print('... working on %s' % file)
        count = np.loadtxt(file, dtype=str)
        indx = np.where([gene_list[i]==count[:,0] for i in range(len(gene_list))])[1]
        count_mtx = np.hstack((count_mtx, count[indx,1].reshape(-1,1)))
        # append sample column heading name 
        sample_header = re.sub(".fastq.gz|.fastq", "", file)
        header.append(sample_header)

    # combine count matrix and column headings
    count_mtx = np.vstack((header, count_mtx))
    
    return count_mtx

def main(argv):
    parsed = parse_args(argv)
    if not os.path.exists(parsed.query):
        sys.exit('ERROR: %s does not exist.' % parsed.input_lookup)
    if not os.path.exists(parsed.gene_list):
        sys.exit('ERROR: %s does not exist.' % parsed.gene_list)
	
    count_matrix = create_count_matrix(parsed.input_lookup, parsed.gene_list)
    np.savetxt(parsed.output_count_file, count_matrix, delimiter=',', fmt='%s')
    print(count_matrix[1:10])	

if __name__ == '__main__':
    main(sys.argv)
