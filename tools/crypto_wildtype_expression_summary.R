#!/usr/bin/env Rscript

# Title     : cyrpto_wildtype_expression_summary.R
# Objective : create .csv with quantile 20 and 99 for evaluation of pertubation experiments
# Created by: chasem@wustl.edu chase.mateusiak@gmail.com
# Created on: 5/8/20

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))

main = function(args){
  # read in cmd line arguments
  parsed = args
  log2_cpm_path = parsed$log2cpm
  quality_summary_path = parsed$quality_assess_2
  query_path = parsed$quality_assess_2
  
  # read in as data.frames
  log2_cpm_df = read.csv(log2_cpm_path, check.names = FALSE) # check.names = false is b/c there are colnames that begin with numbers. without this flag, a X will be added
  query_df_full = read.csv(query_path)
  quality_summary_df_full = readxl::read_excel(quality_summary_path)
  
  # add gene_id to log2_cpm_df
  rownames(log2_cpm_df) = log2_cpm_df$gene_id
  
  # sequentially renumber replicate groups
  replicate_groups = c('GENOTYPE', 'TREATMENT', 'TIMEPOINT')
  query_df_full = splitRenumberRepsCombine(query_df_full, replicate_groups) # arbitrarily (but sequentially within replicate groups) renumber replicates in query_df
  # drop columns that do not have data in query_df_df
  query_df_subset = query_df_full %>%
    select(-LIBRARYPREPARER, -LIBRARYSAMPLENUMBER, -PURPOSE, -S2CDNAPREPARER, -S2CDNASAMPLENUMBER,
           -INDEX1SEQUENCE, -INDEX2NAME, -INDEX2SEQUENCE, -S1CDNAPREPARER, -S1CDNASAMPLENUMBER,
           -RNAPREPARER, -RNASAMPLENUMBER, -HARVESTER, -BIOSAMPLENUMBER, -RNAPREPMETHOD, -ROBOTICRNAPREP,
           -EXPERIMENTDESIGN, -EXPERIMENTOBSERVATIONS, -FLOODMEDIA, -INDUCTIONDELAY)
  # subset quality_summary_df_full
  quality_summary_df_subset = quality_summary_df_full %>%
    select(COUNTFILENAME, STATUS, AUTO_AUDIT, MANUAL_AUDIT, TOTAL, ALIGN_PCT, MUT_FOW)
  # merge quality_summary and query_df to create metadata_df.
  metadata_df = query_df_subset %>%
    left_join(quality_summary_df_subset, by='COUNTFILENAME')%>%
    unite('TREATMENT_TIMEPOINT', c('TREATMENT', 'TIMEPOINT'))
  
  # associate count data with metadata
  expression_meta_df = associateExpressionAndMetadata(log2_cpm_df, metadata_df)
  
  # filter by alignment percentage and total reads, select gene_id, timepoint, treatment
  expression_meta_df = expression_meta_df %>%
    filter(ALIGN_PCT > .79, TOTAL > 1*10^6) %>%
    select(GENE_ID, TREATMENT_TIMEPOINT, LOG2_CPM) %>%
    arrange(GENE_ID, TREATMENT_TIMEPOINT)
  
  wt_log2_cpm_mean_df = expression_meta_df %>%
    group_by(GENE_ID, TREATMENT_TIMEPOINT) %>%
    summarize(quantile = mean(LOG2_CPM))
  
  
} # end main()

parseArguments = function() {
  option_list = list(
    make_option(c('-l', '--log2cpm'),
                help='log2cpm counts created by log2cpm.py'),
    make_option(c('-q', '--standardized_wt_query_sheet'),
                help='result of a search for all wildtypes in metadata database. STANDARDIZED BY DatabaseObject.standardizeDatabaseDataframe()'),
    make_option(c('-s', '--quality_assess_2'),
                help='result of quality_assess_2.py'),
    make_option(c('-o', '--output_FULL_path'),
                help='path to file (full, from current directory through the filename and extension'))
  args = parse_args(OptionParser(option_list=option_list))
  return(args)
} # end parseArguments()

associateExpressionAndMetadata = function(expression_log2_cpm_df, metadata_df){
  #' function to merge
  #'
  #'
  #'
  
  tidy_expression_with_metadata = expression_log2_cpm_df %>%
    gather("COUNTFILENAME", "LOG2_CPM", -gene_id) %>% # reshape data vertically. columns named gene_id, COUNTFILENAME, EXPERSSION vertically
    arrange(gene_id) %>% # sort by gene_id CNAG_00001 lined up in all samples, CNAG_00002, etc
    inner_join(metadata_df, by = "COUNTFILENAME") %>% # inner join with the metadata query. This is inner join, so only common returned
    select_if(~sum(!is.na(.)) > 0) # drop columns with only nas credit: https://stackoverflow.com/a/45383054/9708266
  
  # cast all to upper case
  colnames(tidy_expression_with_metadata) = toupper(colnames(tidy_expression_with_metadata))
  # factor replicate
  tidy_expression_with_metadata$REPLICATE = factor(tidy_expression_with_metadata$REPLICATE)
  
  return(tidy_expression_with_metadata)
} # end associateExpressionAndMetadata()


main(parseArguments()) # call main method