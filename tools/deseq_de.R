#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(DESeq2))
suppressMessages(library(tidyverse))
suppressMessages(library(BiocParallel))
register(MulticoreParam(10))

main = function(parsed_cmd_line_args){
  # main method of script
  
  print('...Parsing cmd line arguments')
  raw_counts_df_path = parsed_cmd_line_args$raw_counts
  metadata_df_path = parsed_cmd_line_args$metadata
  ruvr_unwanted_covariation_path = parsed_cmd_line_args$ruvr_unwanted_covaration_path
  num_unwanted_covariates = as.double(parsed_cmd_line_args$num_unwanted_covariates)
  protein_coding_gene_path = parsed_cmd_line_args$protein_coding_gene_path
  genotype_results_flag = parsed_cmd_line_args$genotype_results_flag
  output_dir = parsed_cmd_line_args$output_directory
  output_name = parsed_cmd_line_args$name
  intercept_only_flag = parsed_cmd_line_args$intercept_only_flag
  # error if no tilde in first position TODO!!
  design_formula = formula(parsed_cmd_line_args$design_formula)
  
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
  metadata_df = factorFormulaColumnsInMetadata(design_formula, metadata_df)
  
  print('...constructing design matrix')
  model_matrix = createModelMatrix(design_formula, metadata_df, intercept_only_flag)
  if (!(is.null(ruvr_unwanted_covariation_path) | is.null(num_unwanted_covariates))){
    
    unwanted_covariate_matrix = as.matrix(read_csv(ruvr_unwanted_covariation_path))
    model_matrix = cbind(model_matrix, unwanted_covariate_matrix[,1:num_unwanted_covariates])
   
  }
  writeOutDataframe(output_path, 'model_matrix', as_tibble(model_matrix))
  
  print('...construct deseq model')
  dds = createDeseqDataObject(raw_counts_df, metadata_df, model_matrix)
  deseq_model = generateDeseqModel(dds)
  deseq_model_path = paste(output_path, 'deseq_model.rds', sep='/')
  saveRDS(deseq_model, deseq_model_path)
  
  
  
} # end main()

convertDesignFormulaToColumnList = function(design_formula){
  
  # extract columns from the design formula as a list
  formula_str = as.character(design_formula)
  column_list = str_split(formula_str, '\\+')
  column_list = lapply(column_list, trimws)[-1]
  
  return(column_list)

} # end convertDesignFormulaToColumnList()

factorFormulaColumnsInMetadata = function(design_formula, df){
  
  # extract columns from the design formula as a list
  column_list = convertDesignFormulaToColumnList(design_formula)
  df[unlist(column_list)] = lapply(df[unlist(column_list)], factor)
  
  if ('LIBRARYPROTOCOL' %in% unlist(column_list)){
    df$LIBRARYPROTOCOL = relevel(df$LIBRARYPROTOCOL, ref='E7420L')
  }
  if ('LIBRARYDATE' %in% unlist(column_list)){
    df$LIBRARYDATE = relevel(df$LIBRARYDATE, ref=max(levels(df$LIBRARYDATE)))
  }
  
  return(df)

} # end factorFormulaColumnsInMetadata

createModelMatrix = function(design_formula, metadata_df, intercept_only_flag){
  
  column_list = convertDesignFormulaToColumnList(design_formula)
  model_matrix = model.matrix(design_formula, metadata_df)
  
  libraryprotocol_librarydate_flag = libraryProtocolDateTest(column_list)
  
  if (libraryprotocol_librarydate_flag){
    col_to_drop = paste0('LIBRARYDATE', min(levels(metadata_df$LIBRARYDATE)))
    col_to_drop_index = match(col_to_drop, colnames(model_matrix))
    model_matrix = model_matrix[,-col_to_drop_index]
  }
  
  if (intercept_only_flag){
    model_matrix = as.matrix(model_matrix[,1])
    colnames(model_matrix) = "Intercept"
  }
  
  return(model_matrix)

} # end createModelMatrix()

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

createDeseqDataObject = function(raw_count_df, metadata_df, model_matrix){
  
  # generate deseq dataset (summarized experiment object) from count, metadata and design formula
  dds = DESeqDataSetFromMatrix(countData = raw_count_df, colData = metadata_df, design = model_matrix)
  
  return(dds)
  
} # end createDeseqDataObject

generateDeseqModel = function(dds){
  
  # construct the deseq model
  deseq_model = DESeq(dds, parallel=TRUE)
  
  return(deseq_model)
  
} # end generateDeseqModel()

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
                help='metadata with all samples corresponding to the columns of the count data x metadata. must include the columns in the design formula'),
    make_option(c('-u', '--ruvr_unwanted_covaration_path'),
                help='path to sheet containing sample x k columns of unwanted variation'),
    make_option(c('-k', '--num_unwanted_covariates'),
                help='number of unwanted covariates to include in the model'),
    make_option(c('-g', '--genotype_results_flag'), action='store_true',
                help='set -g (no input) to write out all genotype results to subdiretory of results directory'),
    make_option(c('-d', '--design_formula'), 
                help='eg ~LIBRARYDATE+GENOTYPE currently does not accept variables with continuous data. base level will be the first level in the factored column(s). Currently not set up for interaction terms'),
    make_option(c('-o', '--output_directory'), 
                help='path to directory to output results'),
    make_option(c('-i', '--intercept_only_flag'), action='store_true', default=FALSE,
                help='set -i to calculate the null model (intercept only). This does require that you put in a valid design formula in -d, but the design formula will not be used.'),
    make_option(c('-n', '--name'),
                help='name of results subdirectory outputed in the path above eg if comparing library date libraryDate_model might be the name'))
  
  args <- parse_args(OptionParser(option_list=option_list))
  return(args)
} # end parseAarguments

main(parseArguments()) # call main method

#for testing
# input_list = list()
# input_list['raw_counts'] = '/home/chase/code/cmatkhan/misc_scripts/deseq_model/data/test_2_counts.csv'
# input_list['metadata'] = '/home/chase/code/cmatkhan/misc_scripts/deseq_model/data/test_2_metadata.csv'
# input_list['design_formula'] = '~LIBRARYDATE'
# input_list['ruvr_unwanted_covaration_path'] = '/home/chase/code/cmatkhan/misc_scripts/deseq_model/results/fullrank_test//unwanted_variation.csv'
# input_list['num_unwanted_covariates'] = 3
# input_list['genotype_results_flag'] = TRUE
# input_list['output_directory'] = '/home/chase/code/cmatkhan/misc_scripts/deseq_model/results'
# input_list['name'] = 'deseq_output_test'

#deviance_df = tibble(gene_id = protein_coding_gene_id_column, deviance_of_fitted_model = mcols(deseq_model)$deviance, saturated_model_deviance = -2*rowSums(dnbinom(counts(deseq_model), mu=counts(deseq_model), size=1/dispersions(deseq_model), log=TRUE)))

# main(input_list)