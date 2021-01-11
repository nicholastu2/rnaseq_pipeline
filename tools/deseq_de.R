#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(DESeq2))
suppressMessages(library(BiocParallel))
register(MulticoreParam(10))

main = function(parsed_cmd_line_args){
  # main method of script
  
  print('...Parsing cmd line arguments')
  dds_path = parsed_cmd_line_args$deseq_data_set
  output_dir = parsed_cmd_line_args$output_directory
  output_name = parsed_cmd_line_args$name
  
  # read in deseq data set as an r object
  dds = readRDS(dds_path)
  
  # run deseq -- this will skip estimating size factors if already provided.
  # it next estimates the dispersion, fits the GLMs and calculates the statistics (wald, default)
  deseq_model = DESeq(dds, parallel=TRUE)
  # save the output of DESeq() as an R object
  deseq_model_path = paste(output_dir, paste0(output_name,'.rds'), sep='/')
  saveRDS(deseq_model, deseq_model_path)
    
} # end main()

parseArguments <- function() {
  # parse and return cmd line input
  
  option_list <- list(
    make_option(c('-d', '--deseq_data_set'),
                help='the dds(deseq data set). see either one of the analysis templates or the deseq docs/vignettes'),
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
# input_list['deseq_data_set'] = '/home/chase/code/cmatkhan/misc_scripts/deseq_model/data/test_2_counts.csv'
# input_list['output_directory'] = '/home/chase/Desktop/tmp/test_results'
# input_list['name'] = 'deseq_output_test'

# main(input_list)
