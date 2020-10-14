library(tidyverse)
library(DESeq2)

# deseq normalize
getDeseqNormCounts = function(count_df, metadata_df){
  metadata_df$lib_prep = as.factor(ifelse(metadata_df$LIBRARYDATE>"2015-10-25", 'stranded', 'unstranded'))
  dds = DESeqDataSetFromMatrix(countData = count_df, colData = metadata_df, design = ~GENOTYPE)
  dds = estimateSizeFactors(dds)
  x = counts(dds, normalized = TRUE)
  return(x)
}

extractResultsLRT = function(counts, metadata, gene_id_column, full_design_matrix, reduced_design_matrix, results_name){
  
  # create the deseq data object
  dds = DESeqDataSetFromMatrix(countData = counts, 
                               colData = metadata, 
                               design=full_design_matrix)
  
  # add genes to the deseq object
  rownames(dds) = gene_id_column
  
  # extract results for effect of timepoint on mutant
  deseq_object = DESeq(dds, test="LRT", reduced = reduced_design_matrix)
  results_object = results(deseq_object, name = results_name)
  
  return(results_object)
  
} # end extractResultsLRT

extractResultsWald = function(counts, metadata, gene_id_column, design_matrix_or_formula, results_name){
  
  if (is.matrix(design_matrix_or_formula)){
    design_matrix_or_formula = design_matrix_or_formula[, colSums(design_matrix_or_formula != 0) > 0]
  }
  
  # create the deseq data object
  dds = DESeqDataSetFromMatrix(countData = counts, 
                               colData = metadata, 
                               design=design_matrix_or_formula,)
  
  # add genes to the deseq object
  rownames(dds) = gene_id_column
  
  # extract results for effect of timepoint on mutant
  deseq_object = DESeq(dds)
  results_object = results(deseq_object, name = results_name)
  
  return(results_object)
  
} # end extractResultsWald

filterCountsByMetadata = function(counts, metadata){
  
  cols_to_keep = colnames(counts) %in% metadata$FASTQFILENAME
  fltr_counts = counts[,cols_to_keep]
  
  return(fltr_counts)
}


