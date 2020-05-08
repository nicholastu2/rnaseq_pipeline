#!/usr/bin/env Rscript

# Title     : log2_counts_per_million
# Objective : convert raw counts to log2 counts per million using edgeR cpm().
#             From the edgeR documentation: If log-values are computed, then a small count, given by prior.count
#             but scaled to be proportional to the library size, is added to y to avoid taking the log of zero.
# Created by: chasem@wustl.edu chase.mateusiak@gmail.com
# Created on: 3/17/20

main = function(args){
  # parse cmd line arguments
  parsed = args
  path_to_raw_counts = parsed$raw_counts
  print('...Creating log2_cpm count matrix...')
  # read in raw counts data.frame
  raw_counts = read.csv(path_to_raw_counts, row.names = 'gene_id', check.names = FALSE) # without check.names = FALSE, R will insert X in front of colnames that start with a number https://stackoverflow.com/a/58951644/9708266
  # convert to edgeR DGEList object
  dgelist = DGEList(raw_counts)
  # cpm returns the log2 of counts per million
  log2_cpm = cpm(dgelist, log=TRUE)
  # convert index row to column 1
  log2_cpm = rownames_to_column(as.data.frame(log2_cpm))
  # rename the converted index row to gene_id
  colnames(log2_cpm)[1] = 'gene_id'

  # write to output path -- NOTE: in cmd line input, the -o is the FULL output path (including filename and extension)
  print('writing log2_cpm matrix to: %s', parsed$output_full_path)
  write_csv(log2_cpm, parsed$output_FULL_path)

} # end main()

suppressMessages(library(optparse))
suppressMessages(library(edgeR))
suppressMessages(library(tidyverse))

parseArguments = function() {
  option_list = list(
    make_option(c('-r', '--raw_counts'),
                help='raw count matrix produced by raw_counts.py'),
    make_option(c('-o', '--output_FULL_path'),
                help='path to file (full, from current directory through the filename and extension'))
  args = parse_args(OptionParser(option_list=option_list))
  return(args)
} # end parseArguments()

main(parseArguments()) # call main method
