#!/usr/bin/env Rscript
suppressMessages(library(optparse))
source('/opt/apps/labs/mblab/software/rnaseq_pipeline/1.0/tools/DE_modules.r')

parse_arguments <- function(){
	option_list <- list(
		make_option(c('-m', '--de_module'), 
					help='Differential analysis module to use. Choose `deseq2` or `edger`'),
		make_option(c('-c', '--count_matrix'), 
					help='Read count matrix (genes x samples).'),
		make_option(c('-d', '--design_table'), 
					help='Design table containing sample grouping indicator.'),
		make_option(c('-q', '--qa_table'), 
					help='QC table containing audit stauts.'),
		make_option(c('-o', '--output_dir'), 
					help='Table of DE genes ranked by ajusted p-value.'))
	opt <- parse_args(OptionParser(option_list=option_list))
	return(opt)
}

## get args and configure paths
parsed_opt <- parse_arguments()
dir.create(parsed_opt$output_dir, showWarnings=FALSE)

## load and parse files
cat('... Parsing metadata\n')
contrast_dict <- parse_metadata(parsed_opt$design_table, parsed_opt$qa_table)
cat('... Loading count matrix\n')
count_matrix <- read.csv(parsed_opt$count_matrix, check.names=FALSE, row.names=1)

## run DE module for each contrast group
if (parsed_opt$de_module == 'deseq2') {
	import_deseq2()
	for (header in names(contrast_dict)) {
		cat('... Working on', header, '\n')
		run_deseq2(count_matrix, contrast_dict, header, parsed_opt$output_dir)
	}
} else if (parsed_opt$de_module == 'edger') {
	import_edger()
	for (header in names(contrast_dict)) {
        cat('... Working on', header, '\n')
        run_edger(count_matrix, contrast_dict, header, parsed_opt$output_dir)
	}
} else {
	cat('ERROR: DE module', parsed_opt$de_module, 'is unavailable.\n')
}

# output DESeq2 QA
#### create heatmap to examine similiarities btwn samples
# rld <- rlogTransformation(ddsTxi, blind=TRUE)
# 
# distsRL <- dist(t(assay(rld)))
# distsRL
# mat <- as.matrix(distsRL)
# rownames(mat) <- colData(rld)$condition
# colnames(mat) <- colData(rld)$condition
# 
# hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(255)
# pheatmap(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
# 
# # PCA
# # from the vignette re: vst : DESeq2 offers transformations for count data that stabilize the variance across the mean: the regularize logarithm (rlog) and the variance stabilizing transformation (VST). 
# # These have slightly different implementations, discussed a bit in the DESeq2 paper and in the vignette, but a similar goal of stablizing the variance across the range of values. 
# # Both produce log2-like values for high counts
# vsd <- vst(ddsTxi)
# plotPCA(vsd, 'condition')
# 
# ### plot dispersion
# ddsTxi_disp <- estimateSizeFactors(ddsTxi)
# ddsTxi_disp <- estimateDispersions(ddsTxi_disp)
# plotDispEsts(ddsTxi_disp, main = "Dispersion")
