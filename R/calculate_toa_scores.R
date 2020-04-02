#' Calculate a tissue-of-action (TOA) scores for each genetic signal
#'
#' This function takes the tissue scores for each individual SNP outputted from
#' calculate_tissue_vectors function and calculates a TOA score for each tissue/cell type for each
#' genetic signal. A signal would correspond to a credible set from genetic fine-mapping or simply
#' an index SNP from GWAS
#'
#' @param snp.tissvec.df Dataframe of tissue scores for each SNP (output from calculate_tissue_vectors function)
#' @export
calculate_toa_scores <- function(snp.tissvec.df){
  '%>%' <- magrittr::'%>%'
  out.df <- c()
  pb <- txtProgressBar(min=0,max=length(unique(snp.tissvec.df$SIGNAL)),style=3)
  count <- 0
  for (sig in unique(snp.tissvec.df$SIGNAL)){
    count <- count+1
    sig.df <- dplyr::filter(snp.tissvec.df,SIGNAL==sig)
    cumsum <- sig.df$VALUE %>% sum(.)
    tiss.matrix <- dplyr::select(sig.df,-one_of("SIGNAL","SNPID","CHR","POS","VALUE")) %>% as.matrix(.)
    accounted.value <- colSums(tiss.matrix) %>% sum(.)
    unclassified.score <- (cumsum - accounted.value) %>% round(.,digits=4)
    sum.df <- colSums(tiss.matrix) %>% round(.,digits=4) %>% t(.) %>% as.data.frame(.)
    sum.df$unclassified <- unclassified.score
    build.df <- data.frame("SIGNAL"=sig,stringsAsFactors = F) %>% cbind(.,sum.df)
    out.df <- rbind(out.df,build.df)
    setTxtProgressBar(pb,count)
  }
  return(out.df)
}
