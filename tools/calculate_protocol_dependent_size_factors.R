#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(DESeq2))
suppressMessages(library(tidyverse))

main = function(parsed_cmd_line_args){
  # main method of script
  
  print('...Parsing cmd line arguments')
  raw_counts_df_path = parsed_cmd_line_args$raw_counts
  metadata_df_path = parsed_cmd_line_args$metadata
  factor_column_list = strsplit(parsed_cmd_line_args$factor_column_list, ",")
  output_dir = parsed_cmd_line_args$output_directory
  output_name = parsed_cmd_line_args$name
  
  # create output directory
  output_path = paste(output_dir,output_name, sep='/')
  print(paste0('output directory created at: ', output_path))
  dir.create(output_path)
  
  print('...reading in raw counts')
  raw_counts_df = read_csv(raw_counts_df_path)
  
  print('...reading in metdata')
  metadata_df = read_csv(metadata_df_path)
  metadata_df$LIBRARYDATE = as.Date(metadata_df$LIBRARYDATE, format="%m.%d.%y")
  
  print('...factoring design formula columns')
  metadata_df = factorFormulaColumnsInMetadata(factor_column_list, metadata_df)
  
  protocol_size_factors = calculateProtocolSizeFactors(metadata_df, raw_counts_df)
  df = as.data.frame(protocol_size_factors)
  rownames(df) = names(protocol_size_factors)
  write.csv(df, paste(output_path, 'protocol_specific_size_factors.csv', sep='/'))
  
  raw_counts_df = orderCountColumns(raw_counts_df, protocol_size_factors)
  write_csv(raw_counts_df, paste(output_path, 'raw_counts_ordered_by_size_factors.csv', sep='/'))
  
  metadata_df = orderMetadataRows(metadata_df, protocol_size_factors)
  write_csv(raw_counts_df, paste(output_path, 'metadata_ordered_by_size_factors.csv', sep='/'))

} # end main()

## NOTE: THE COMMENTED OUT FUNCTION IS THE "NEW" METHOD OF DOING THE DESIGN MATRIX. THE FUNCTION WITH THE NEW
##       PREP IN THE INTERCEPT IS THE "OLD" 2020121.
##       TWO VERSIONS ARE CURRENTLY INCLUDED FOR TESTING
# factorFormulaColumnsInMetadata = function(column_list, df){
#   
#   df[unlist(column_list)] = lapply(df[unlist(column_list)], factor)
#   
#   if ('LIBRARYPROTOCOL' %in% unlist(column_list)){
#     df$LIBRARYPROTOCOL = relevel(df$LIBRARYPROTOCOL, ref='SolexaPrep')
#   }
#   if ('LIBRARYDATE' %in% unlist(column_list)){
#     df$LIBRARYDATE = relevel(df$LIBRARYDATE, ref=min(levels(df$LIBRARYDATE)))
#   }
#   
#   return(df)
# 
# } # end factorFormulaColumnsInMetadata

factorFormulaColumnsInMetadata = function(column_list, df){
  
  df[unlist(column_list)] = lapply(df[unlist(column_list)], factor)
  
  if ('LIBRARYPROTOCOL' %in% unlist(column_list)){
    df$LIBRARYPROTOCOL = relevel(df$LIBRARYPROTOCOL, ref='E7420L')
  }
  if ('LIBRARYDATE' %in% unlist(column_list)){
    df$LIBRARYDATE = relevel(df$LIBRARYDATE, ref=max(levels(df$LIBRARYDATE)))
  }
  
  return(df)
  
} # end factorFormulaColumnsInMetadata

libraryProtocolDateTest = function(column_list){
  
  library_date_flag = FALSE
  if ('LIBRARYDATE' %in% unlist(column_list)) {
    library_date_flag = TRUE
  }
  
  library_protocol_date_flag = FALSE
  if (library_date_flag == TRUE){
    if ('LIBRARYPROTOCOL' %in% unlist(column_list)){
      library_protocol_date_flag = TRUE
    }
  }
  
  return(library_protocol_date_flag)
  
} # end libraryProtocolDateTest()

calculateProtocolSizeFactors = function(metadata_df, raw_counts){
  #' new function 20201123 to calculate protocol specific size factors -- needs editing/debugging/error handling in script
  #' scripting here is copied from notebook -- needs cleaning up
  
  old_library_metadata_df = metadata_df %>% filter(LIBRARYPROTOCOL == "SolexaPrep")
  new_library_metadata_df = metadata_df %>% filter(LIBRARYPROTOCOL != "SolexaPrep")
  
  old_library_counts = as_tibble(raw_counts) %>% select(old_library_metadata_df$FASTQFILENAME)
  new_library_counts = as_tibble(raw_counts) %>% select(new_library_metadata_df$FASTQFILENAME)
  
  dds_old = DESeqDataSetFromMatrix(colData=old_library_metadata_df, countData=old_library_counts, design=~1)
  dds_old = estimateSizeFactors(dds_old)
  old_size_factors = sizeFactors(dds_old)
  
  dds_new = DESeqDataSetFromMatrix(colData=new_library_metadata_df, countData=new_library_counts, design=~1)
  dds_new = estimateSizeFactors(dds_new)
  new_size_factors = sizeFactors(dds_new)
  
  size_factors_by_protocol = c(old_size_factors, new_size_factors)
  
  return(size_factors_by_protocol)
  
} # end calculateProtocolSizeFactors()

orderCountColumns = function(raw_counts, size_factors_by_protocol){
  #' order columns of counts based on order of size factors
  counts_ordered = raw_counts[names(size_factors_by_protocol)]
  counts_ordered = replace(counts_ordered, is.na(counts_ordered), 0)
  
  return(counts_ordered)
  
} # end orderCountColumns()

orderMetadataRows = function(metadata_df, size_factors_by_protocol){
  #' order rows of metadata based on order of size factors
  
  metadata_df_ordered = metadata_df[match(names(size_factors_by_protocol), metadata_df$FASTQFILENAME),]
  
  return(metadata_df_ordered)
  
} # end orderMetadataRows()

writeOutDataframe = function(output_path, chart_name, df){
  
  # output path
  csv_output_path = paste(output_path, paste0(chart_name, '.csv'), sep='/')
  # tell user whats what  
  print(paste0('writing sheet: ', csv_output_path))
  # write
  write_csv(df, csv_output_path)
  
} # end writeOutDataframe()

parseArguments <- function() {
  # parse and return cmd line input
  
  option_list <- list(
    make_option(c('-r', '--raw_counts'),
                help='raw count matrix (genes x samples)'),
    make_option(c('-m', '--metadata'), 
                help='metadata with all samples corresponding to the columns of the count data x metadata. 
                      must include the columns in the design formula'),
    make_option(c('-f', '--factor_column_list'), 
                help='comma separated list NO SPACES of columns to factor, eg LIBRARYDATE, GENOTYPE'),
    make_option(c('-o', '--output_directory'), 
                help='path to directory to output results'),
    make_option(c('-n', '--name'),
                help='name of results subdirectory outputed in the path above eg if comparing library date libraryDate_model might be the name'))
  
  args <- parse_args(OptionParser(option_list=option_list))
  return(args)
} # end parseAarguments

main(parseArguments()) # call main method

# for testing
# input_list = list()
# input_list['raw_counts'] = '/home/chase/code/cmatkhan/misc_scripts/deseq_model/data/test_2_counts.csv'
# input_list['metadata'] = '/home/chase/code/cmatkhan/misc_scripts/deseq_model/data/test_2_metadata.csv'
# input_list['size_factors'] = '/mnt/htcf_scratch/chasem/rnaseq_pipeline/experiments/size_factors/data/size_factors.csv'
# input_list['factor_column_list'] = 'LIBRARYDATE,GENOTYPE'
# input_list['design_matrix'] = '/home/chase/code/cmatkhan/misc_scripts/deseq_model/data/librarydate_genotype_model_matrix.csv'
# input_list['genotype_results_flag'] = TRUE
# input_list['output_directory'] = '/home/chase/Desktop/tmp/test_results'
# input_list['name'] = 'deseq_output_test'

#deviance_df = tibble(gene_id = protein_coding_gene_id_column, deviance_of_fitted_model = mcols(deseq_model)$deviance, saturated_model_deviance = -2*rowSums(dnbinom(counts(deseq_model), mu=counts(deseq_model), size=1/dispersions(deseq_model), log=TRUE)))

# main(input_list)