library(tidyverse)
library(DESeq2)

readInDeseqModel = function(path_to_deseq_model_rds, gene_ids){
  #' the gene_id column should be unnecessary b/c it should be included in the calculation of the deseq_model object
  #' However, until that is the case, must be provided
  
  dds = readRDS(path_to_deseq_model_rds)
  rownames(dds) = gene_ids
  
  return(dds)
  
} # end readInDeseqModel()

predictedModelMeans = function(deseq_model){
  #' This method extracts the predicted means of each gene in each sample. See the 2014 deseq paper, equations
  #' 1 and 2 in the methods section.
  #' s is the size factors, by default calculated per sample, calculated by deseq (see paper). extract from deseq_model
  #' with sizeFactors(deseq_model)
  #' variable are named according to the methods section, equations 1 and 2, of the 2014 deseq paper
  #' NOTE: returns mu in scale of raw_counts
  
  # model matrix is m=sample  n=model predictors, coefficient matrix is m=gene n=predictors
  log_q = coef(deseq_model) %*% t(design(deseq_model))
  # q is i=gene j=samples
  q = 2^log_q # exponentiate to un-do the log2 transformation
  mu = t(t(q)*sizeFactors(deseq_model)) # these are now the predicted counts in the scale of the raw_counts
  
  # to test the above multiplication
  # mdat <- matrix(c(1,2,3, 11,12,13), nrow = 3, ncol = 2, byrow=FALSE)
  #t(t(mdat) * seq(1,2))
  
  # name rows and columns
  rownames(mu) = rownames(deseq_model)
  colnames(mu) = colnames(counts(deseq_model))
  
  return(mu)
  
} # end calculateModelPredictions()

modelCountLogProbability = function(deseq_model, mu_matrix){
  #' calculate the likelihood of K_ij given NB(mu_ij, dispersion_i)
  #' see 2014 DESeq2 paper, equations 1 and 2
  #' probabilities are log10 transformed
  
  # return probability of observing K_ij given NB(mu_ij, dispersion_i)  
  log_probability_matrix = dnbinom(counts(deseq_model, normalized=FALSE), # normalized=FALSE is default, included to be explicit
                                   mu=mu_matrix, 
                                   size=1/dispersions(deseq_model), # this formulation of dnbinom is adapted from Michael Love here: http://supportupgrade.bioconductor.org/p/117448/ 
                                   log=TRUE) # natural logarithm
  
  return(log_probability_matrix)
} # end modelCountLogProbability()

nullModelMeanMatrix = function(deseq_model){
  #' for purposes of calculating null model using modelCountLogProbability
  #' create matrix that is the rowMeans (genewise means) repeated by
  #' the number of columns in the count matrix
  
  counts = counts(deseq_model, normalized=TRUE)
  row_means = rowMeans(counts)
  
  
} # end nullModelMeanMatrix()

proposedModelNegTwoLogLikelihood = function(deseq_model){
  
  saturated_model_neg_loge_likelihood = -2*rowSums(modelCountLogProbability(deseq_model, predictedModelMeans(deseq_model)))
  
  return(saturated_model_neg_loge_likelihood)
}

saturatedModelNegTwoLogLiklihood = function(deseq_model){
  #'  The below function is based on a response from Michael love here:
  #'      http://supportupgrade.bioconductor.org/p/117448/   
  
  saturated_model_neg_loge_likelihood = -2*rowSums(modelCountLogProbability(deseq_model, counts(deseq_model)))
  
  return(saturated_model_neg_loge_likelihood)
  
}# end modelCountLogProbability()

calculateDevianceResiduals = function(deseq_model){
  #' NOTE: deviance residuals = 
  #' while Residual Deviance of a model = 
  proposed_count_probability_df = modelCountLogProbability(deseq_model, predictedModelMeans(deseq_model))
  
  saturated_model_probability_df = modelCountLogProbability(deseq_model, counts(deseq_model))
  
  # deviance residuals of each sample are 2(log(p(k_ij|proposed)) - log(p(k_ij|saturated))
  diff_in_log_probability_df = 2* (saturated_model_probability_df - proposed_count_probability_df)
  # take the square root of the exponentiated elements
  unsigned_deviance_residuals = sqrt(diff_in_log_probability_df) # in log space
  # determine sign based on K_ij - q_ij
  sign_df = sign(counts(deseq_model, normalized=FALSE) - predictedModelMeans(deseq_model))
  # multiply the sign_df by the unsigned_deviance_residuals NOTE: NOT MATRIX MULTIPLICATION
  signed_deviance_residuals = unsigned_deviance_residuals * sign_df
  
  return(signed_deviance_residuals) # in log space
  
}

#calculateDevianceResiduals(deseq_librarydate)

calculateResidualDeviance = function(deseq_model){
  #' Deviance Residuals are:
  #'     Residual Deviance = 2(log_likelihood(Saturated_Model) - log_likelihood(Proposed_Model))
  #'           (this is analogous to residuals in ordinary least squares per statsquest)
  #'  The below function is based on a response from Michael love here:
  #'      http://supportupgrade.bioconductor.org/p/117448/   
  
  # per michael love: The deviance as calculated by DESeq2 is -2 times the log likelihood
  deviance_of_fitted_model = mcols(deseq_model)$deviance
  
  deviance_of_saturated_model = saturatedModelNegTwoLogLiklihood(deseq_model)
  # note: each model multiplied by -2. hence switch in order
  residual_deviance = deviance_of_fitted_model - deviance_of_saturated_model
  
  return(residual_deviance)
} # end calculateResidualDeviance()

calculateVarianceCovariance = function(deviance_residuals_df){
  
  return(var(deviance_residuals_df, na.rm=TRUE))
}