
plotNatCoverageVsExpression = function(meta_qual_df, log2_cpm_df, title_prefix){
  
  log2_meta_qual = 
  
  # plot coverage x log2cpm for NAT, color by whether or not NAT is expected
  g = meta_qual_df %>%
    ggplot(aes(NAT_COVERAGE, NAT_LOG2CPM, color=MARKER_1=='NAT', shape=AUTO_AUDIT))+
    geom_point(alpha=.5, size=3)+
    ggtitle(paste0(title_prefix, 'coverage vs log2cpm for NAT, colored by NAT expected'))
  
  return(g)

 }  
  
  
plotNatCoverageVsExpression = function(meta_qual_df, log_cpm_df, title_prefix){
  
  g = meta_qual_df %>%
    ggplot(aes(G418_COVERAGE, G418_LOG2CPM, color=MARKER_1=='G418', shape=AUTO_AUDIT))+
    geom_point(alpha=.5, size=3)+
    ggtitle(paste0(title_prefix, 'coverage vs log2cpm for G418, colored by G418 expected'))
  
  return(g)
  
}  
  
  
  
  
