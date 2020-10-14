library(tidyverse)
library(DESeq2)

createFullRLETable = function(counts_df, gene_id_column){
  counts_df = as_tibble(counts_df)
  counts_df$gene_id = gene_id_column
  
  # log2 the table (add pseudo count of 1)
  log2_counts_df = log2(counts_df %>% select(-gene_id) + 1)
  
  # calculate median expression for each gene across samples
  all_gene_medians = apply(log2_counts_df, 1, median)
  
  # calculate deviations from the median
  rle_table_full = sweep(log2_counts_df, 1, all_gene_medians, '-')
  
  return(rle_table_full)
  
} # end createFullRLETable()

rleSummary = function(rle_table_full){
  # calculate median deviation by sample
  median_deviation_by_sample = apply(rle_table_full, 2, median)
  
  # calculate twenty fifth percentile of deviations from median
  q1 = apply(rle_table_full, 2, quantile, .25)
  # ditto seventy fifth
  q3 = apply(rle_table_full, 2, quantile, .75)
  # calculate interquartile range
  iqr = q3 - q1
  # min/max whisker threshold = q1/3 +/- 1.5*iqr
  min_whisker_threshold = q1 - 1.5*iqr
  
  min_whisker_value = c()
  for (i in seq(1,length(colnames(rle_table_full)))){
    x = min(rle_table_full[,i][rle_table_full[,i]> min_whisker_threshold[i]])
    min_whisker_value = c(min_whisker_value,x)
  }
  
  max_whisker_threshold = q3 + 1.5*iqr
  max_whisker_value = c()
  for (i in seq(1,length(colnames(rle_table_full)))){
    x = max(rle_table_full[,i][rle_table_full[,i]< max_whisker_threshold[i]])
    max_whisker_value = c(max_whisker_value,x)
  }
  
  inter_whisker_range = max_whisker_value - min_whisker_value
  
  num_outliers_below =c()
  for (i in seq(1,length(colnames(rle_table_full)))){
    x = length(rle_table_full[,i][rle_table_full[,i] < min_whisker_threshold[i]])
    num_outliers_below = c(num_outliers_below,x)
  }
  
  most_extreme_below = apply(rle_table_full,2,min)
  for (i in seq(1,length(most_extreme_below))){
    if (most_extreme_below[i] >= min_whisker_threshold[i]){
      most_extreme_below[i] = NA
    }
  }
  
  num_outliers_above =c()
  for (i in seq(1,length(colnames(rle_table_full)))){
    x = length(rle_table_full[,i][rle_table_full[,i] > max_whisker_threshold[i]])
    num_outliers_above = c(num_outliers_above,x)
  }
  
  most_extreme_above = apply(rle_table_full,2,max)
  for (i in seq(1,length(most_extreme_above))){
    if (most_extreme_above[i] <= max_whisker_threshold[i]){
      most_extreme_above[i] = NA
    }
  }
  
  rle_table_summary = tibble(FASTQFILENAME = colnames(rle_table_full), SAMPLE_DEVIATION_MEDIAN = median_deviation_by_sample, ABS_SAMPLE_DEVIATION_MEDIAN = abs(median_deviation_by_sample), TWENTY_FIVE_QUANTILE = q1, SEVENTY_FIVE_QUANTILE = q3, INTERQUARTILE_RANGE = iqr, MIN_WHISKER_VALUE = min_whisker_value, NUM_OUTLIERS_BELOW = num_outliers_below, MOST_EXTREME_BELOW = most_extreme_below, MAX_WHISKER_VALUE = max_whisker_value, NUM_OUTLIER_ABOVE = num_outliers_above, MOST_EXTREME_ABOVE = most_extreme_above)
  
  return (rle_table_summary)
  
} # end rleSummary()

rleBoxplots= function(rle_table_full, meta_qual_df, table_name, feature){
  #' rle_table_full output of rleTableFull() above
  #' meta_qual_df is a metadata + quality assess dataframe
  #' table_name is the name of the graph
  #' part of the table name, the feature 
  
  df = stack(rle_table_full)
  
  stacked_rle_meta_qual_df = df %>% left_join(meta_qual_df, by=c('ind'='FASTQFILENAME'))
  
  ylim1 = boxplot.stats(df$values)$stats[c(1, 5)]
  
  # order by different variables in metadata (eg library date) <-- TODO: make this option
  
  # rle_boxplots = ggplot(stacked_rle_meta_qual_df, aes(reorder(ind, PROTEIN_CODING_COUNTED), values))+
  #   geom_boxplot(outlier.shape=NA)+
  #   coord_cartesian(ylim = ylim1*2.5)+
  #   theme(axis.text.x=element_blank(),
  #         axis.ticks.x=element_blank())+
  #   xlab('Samples')+
  #   ylab('Deviation from Median')+
  #   ggtitle('90 Minute Induction RLE plots')
  # 
  # plot(rle_boxplots)
  # 
  
  rle_boxplots = ggplot(stacked_rle_meta_qual_df, aes(ind, values))+
    geom_boxplot(outlier.shape=NA, aes(fill=GENOTYPE))+
    coord_cartesian(ylim = ylim1*2.5)+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    xlab('Samples')+
    ylab('Deviation from Median')+
    geom_hline(yintercept = c(.7,-.7))+
    ggtitle(paste0(table_name, '_', feature, ' 90 Minute Induction RLE plots'))
  
  return(rle_boxplots)
  
} # end rleBoxPlots()