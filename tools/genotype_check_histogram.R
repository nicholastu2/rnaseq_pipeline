#!/usr/bin/env Rscript

# Title     : genotype_check_histogram
# Objective : create a histogram(s) of  median(log2cpm) of replicate set
#             with perturbed gene(s) from each replicate and drug markers overlaid. See templates for example.
# Created by: chase mateusiak chasem@wustl.edu chase.mateusiak@gmail.com
# Created on: 3/17/20

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
#suppressMessages(library(calecopal)) # this is a package of color palettes that can be found here:


main = function(){
  parsed = parseArguments()
  # parse cmd line arguments
  parsed = parseArguments()
  path_to_query_sheet = parsed$standardized_query_sheet
  path_to_log2_cpm = parsed$log2_cpm
  input_cols = toupper(str_split(parsed$experiment_columns, ',')[[1]]) # parse on comma and cast to uppercase
  output_dir = parsed$output

  print(path_to_query_sheet)
  print(path_to_log2_cpm)
  print(input_cols)
  print(output_dir)

  # # read in data
  # log2_cpm = read_csv(path_to_log2_cpm)
  # # read in query and raw counts
  # parsed_query_df = read_csv(path_to_query_sheet)
  #
  # # create list of columns to parse out of query_sheet
  # columns_to_parse = c('COUNTFILENAME', 'REPLICATE', input_cols)
  #
  # # split the query_df into separate sheets based on cols_of_interest
  # split_parsed_query_df = splitParsedQuery(parsed_query_df, columns_to_parse, input_cols)
  #
  # createHistograms(split_parsed_query_df, input_cols, log2_cpm, output_dir)

}  # end main()

parseArguments = function() {
	option_list = list(
		make_option(c('-q', '--standardized_query_sheet'),
					help='a query sheet for the samples that has been standardized by class StandardDataFormat'),
        make_option(c('-c', '--log2_cpm'),
					help='raw count matrix produced by raw_counts.py'),
        make_option(c('-e', '--experiment_columns'),
					help='columns from the query sheet to be used to separate replicate groups. Eg for crypto, -col genotype treatment timepoint'),
		make_option(c('-o', '--output'),
					help='path to directory in which to deposit graphs. No trailing /'))
	args = parse_args(OptionParser(option_list=option_list))
	return(args)
} # end parseArguments()

splitParsedQuery = function(parsed_query_df, columns_to_parse, input_cols){

  # split parsed_query_df into groups based on cols_of_interest
  split_parsed_query_df = parsed_query_df %>%
                          select(columns_to_parse) %>%
                          # The ugly code below 'unquoting' the character vector in order to pass
                          # the arguments as column variables see https://adv-r.hadley.nz/quasiquotation.html
                          group_split(!!! rlang::syms(input_cols))

  return(split_parsed_query_df)

} # end splitParsedQuery

createHistograms = function(split_parsed_query_df, input_cols, log2_count_cpm, output_dir){

  for(group in split_parsed_query_df){

    group_desc = paste(as.character(unique(group[,input_cols])), collapse = "_")
    x_label = "median of log2(cpm) of replicates"
    y_label = "frequency"

    # arbitrarily re-number replicates
    group$REPLICATE = 1:nrow(group) # this should be done already in StandardizeDataFormat

    # subset log2_count_cpm by the countfilename (samples) in a given group in split_parsed_query_df
    subset_log2_cpm = log2_count_cpm[, group$COUNTFILENAME]
    subset_log2_cpm = as.data.frame(subset_log2_cpm)
    subset_log2_cpm$MEDIAN = apply(subset_log2_cpm, 1, median)
    # summarize by median
    subset_log2_cpm = select(subset_log2_cpm, median)

    # plot histogram of median log2 cpm of a given group in split_parsed_query_df
    overall_dist = ggplot(subset_log2_cpm, aes(median))+
      geom_histogram(na.rm = TRUE, bins = 105, alpha = .5, fill = 'tan2')

    # extract genotype of the group (there will be only one)
    genotype = unique(group[,'GENOTYPE'])[[1]]
    # remove "_over" if it is present
    if(str_detect(genotype, "_over")){
      genotype = str_remove(genotype, "_over")
    }

    countfilename_list = as.character(group$COUNTFILENAME)
    if(str_detect(genotype, "\\.")){
      genotype = str_split(genotype, "\\.")[[1]]
    }

    # add drug markers
    drug_markers = c("CNAG_G418", "CNAG_NAT")

    # extract expression of perturbed gene
    perturbed_df = as_tibble(log2_count_cpm[genotype, countfilename_list])
    if (ncol(perturbed_df) == 1){
      perturbed_df$SAMPLE = countfilename_list
      perturbed_df$GENOTYPE = genotype
    } else {
      perturbed_df = stack(perturbed_df)
      perturbed_df$GENE = genotype
      colnames(perturbed_df)[2] = "SAMPLE"
      colnames(perturbed_df)[3] = "GENOTYPE"
    }
    colnames(perturbed_df) = toupper(colnames(perturbed_df))

    # get expression of drug markers
    drug_markers = c("CNAG_G418", "CNAG_NAT")
    drug_markers_df = as_tibble(log2_count_cpm[drug_markers, countfilename_list])
    drug_markers_df = stack(drug_markers_df)
    drug_markers_df$gene = drug_markers
    colnames(drug_markers_df)[2] = "SAMPLE"
    colnames(drug_markers_df)[3] = "GENOTYPE"
    colnames(drug_markers_df) = toupper(drug_markers_df)

    # plot points representing expression of a genotype (or marker) in a given sample
    overall_dist = overall_dist +
      # uncomment this to display dots by sample (countfilename)
      #geom_point(data = perturbed_df, aes(perturbed_df$value,y=0, color = perturbed_df$sample)) +
      geom_text_repel(data = perturbed_df, aes(x = perturbed_df$VALUE, y = perturbed_df$VALUE, label = perturbed_df$GENOTYPE),
                      nudge_y = 200,
                      vjust = 0,
                      segment.size = .2) +
      geom_point(data = drug_markers_df, aes(drug_markers_df$VALUE,y=100, color = drug_markers_df$GENOTYPE))+
      theme(legend.justification = c(1,1), legend.position = c(1,1))+
      labs(title=group_desc, y=y_label, x=x_label)+
      xlim(0,15)+
      ylim(0,400)

    output_path = file.path(output_dir, paste(group_desc, '.png', sep=''))
    ggsave(output_path, plot=overall_dist)
  } # end for
} # end createHistograms

main() # call main method