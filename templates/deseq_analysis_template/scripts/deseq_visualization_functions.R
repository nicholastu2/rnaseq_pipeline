library(tidyverse)
library(DESeq2)
library(factoextra)
library(ggcorrplot)

# to install ggcorrplot
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggcorrplot")
#http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2


statistics_functions_script = '~/projects/brentlab/deseq_model/deseq_statistics_functions.R'
source(statistics_functions_script)

plotResidualDeviance = function(deseq_model){
  #' plot residual deviance by gene (gene on x, residual deviance of the gene on y)
  #' and a histogram of residual deviance
  
  residual_deviance = calculateResidualDeviance(deseq_model)
  plot(residual_deviance)
  hist(residual_deviance)
} # end plotResidualDeviance()

calculatePrincipalComponents = function(residuals_df, metadata_df, scree_title){
  #' calculate principle components and plot variance explained (scree plot)
  #' :params residuals_df: residuals of the deseq model, preferrably logged
  #' :params metadata_df: the metadata used for the model
  #' :params scree_title: currently useless. Does not actually appear on plot.
  #' :returns: the princple components of the residuals joined with the metadata
  #' :side effects: plots the scree plot
  
  # calculate principal components on residuals after log2 transformation (add psueocount)
  no_na_resid_df = residuals_df[is.na(residuals_df)] = 0
  no_na_resid_df = residuals_df[is.infinite(residuals_df)] = 0
  residuals_prcomp_object = prcomp(residuals_df)
  # scree plot
  plot(fviz_eig(residuals_prcomp_object), main = scree_title)
  
  # convert to tibble
  residuals_pc_df = as_tibble(residuals_prcomp_object$rotation)
  # add rownames
  residuals_pc_df$FASTQFILENAME = rownames(residuals_prcomp_object$rotation)
  # join with metadata
  residuals_pc_df = dplyr::inner_join(residuals_pc_df, metadata_df, on=FASTQFILENAME)
  
  return(residuals_pc_df)
  
} # end calculatePrincipalComponents

residualPcaPlot = function(log_space_residual_pc_df, graph_title){
  #' plot pca of principle components of residuals (preferrably in log space)
  #' :params log_space_residual_pd_df: the (preferrably logged) pc of the residuals joined with metadata
  #' :params graph_title: a title for the graph (eg, the design formula for the model that produced the residuals)
  #' :returns: None. plots the charts

  # TODO: GENERALIZE TO CREATE LIST OF PLOTS -- PASS LIST OF COLUMN VARIABLES IN TO PLOT
  # note: design_formula passed as a formula object
  
  # factor library_date
  log_space_residual_pc_df$LIBRARYDATE = as.factor(log_space_residual_pc_df$LIBRARYDATE)
  
  g_librarydate = ggplot(log_space_residual_pc_df, aes(PC1,PC2))+geom_point(aes(color=LIBRARYDATE))+ggtitle(graph_title)
  plot(g_librarydate)
  g_libraryprotocol = ggplot(log_space_residual_pc_df, aes(PC1,PC2))+geom_point(aes(color=LIBRARYPROTOCOL))+ggtitle(graph_title)
  plot(g_libraryprotocol)
  g_genotype = ggplot(log_space_residual_pc_df, aes(PC1,PC2))+geom_point(aes(color=GENOTYPE))+ggtitle(graph_title)
  plot(g_genotype)
  g_timepoint = ggplot(log_space_residual_pc_df, aes(PC1,PC2))+geom_point(aes(color=TIMEPOINT), size=4)+
    ggtitle(graph_title)
  plot(g_timepoint)
  
  g_combined = ggplot(log_space_residual_pc_df, aes(PC1,PC2))+geom_point(aes(color=LIBRARYDATE, shape=TIMEPOINT), size=4)+
    geom_point(data=log_space_residual_pc_df[log_space_residual_pc_df$GENOTYPE == 'CNAG_00000',], color="black", size=1.5)+
    ggtitle(paste(graph_title, 'black dots == wt', sep=' '))
  plot(g_combined)
  
  # library_date_output = paste(output_path, 'residual_pca_by_library_date.pdf', sep='/')
  # print('writing library_date PCA')
  # ggsave(filename = library_date_output, plot = g_librarydate, device='pdf', height=8, width=12)
  # library_protocol_output = paste(output_path, 'residual_pca_by_library_protocol.pdf', sep='/')
  # print('writing library_protocol PCA')
  # ggsave(filename = library_protocol_output, plot = g_libraryprotocol, device='pdf', height=8, width=12)
  # genotype_output = paste(output_path, 'residual_pca_by_genotype.pdf', sep='/')
  # print('writing genotype PCA')
  # ggsave(filename = genotype_output, plot = g_genotype, device='pdf', height=8, width=12)
  
} # end createPcaPlots()

plotVarianceByModel = function(){
  #' this is a todo: write a function to plot variance exlpained by nested models on same scale
  
}

pValueHistograms = function(res, title){
  #' given a deseq results table and a title, plot a histogram of the p-values
  #' Code is copied from the deseq2 vignette
  #' :params res: a deseq object results table
  #' :params title: a title for the graph
  #' :returns: None. Plots out the chart
  
  use <- res$baseMean > metadata(res)$filterThreshold
  h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
  h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
  colori <- c(`do not pass`="khaki", `pass`="powderblue")
  barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
          col = colori, space = 0, main = title, ylab="frequency")
  text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
       adj = c(0.5,1.7), xpd=NA)
  legend("topright", fill=rev(colori), legend=rev(names(colori)))
}