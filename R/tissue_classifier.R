#' Classify genetic signals using a simple rule-based classifier
#'
#' This function takes the calculated toa scores for each genetic signal and applies a simple
#' rule-based classifier to assign each signal to a tissue/cell-type. The signal is assigned to the
#' tissue/cell-type that has the highest TOA score and exceeds a user-specified score threshold (default 0.2).
#' The user can also specify a threshold for determining shared signals. If two or more tissues yield TOA scores
#' that fall within this specified range of each other, the signal is designated as a "shared" signal. If
#' designated as shared, the responsible tissues will be indicated in the output data frame.
#'
#' @param toa.df Dataframe of tissue scores for each SNP (output from calculate_tissue_vectors function)
#' @param tissue_threshold Dataframe of tissue scores for each SNP (output from calculate_tissue_vectors function)
#' @param shared_threshold Dataframe of tissue scores for each SNP (output from calculate_tissue_vectors function)
#' @export
tissue_classifer <- function(toa.df,tissue_threshold=0.2,shared_threshold=0.1){
  pb <- txtProgressBar(min=1,max=dim(toa.df)[1],style=3)
  out.df <- c()
  for (i in 1:dim(toa.df)[1]){
    row.df <- toa.df[i,]
    tiss.scores <- row.df %>% dplyr::select(.,-one_of("SIGNAL","unclassified")) %>%
      sort(.,decreasing = TRUE)
    tiss.names <- names(tiss.scores)
    keep.scores <- tiss.scores[tiss.scores > tissue_threshold]
    keep.names <- tiss.names[tiss.scores > tissue_threshold]
    if (length(keep.scores)==0){
      classification <- "unclassified"
      tissues <- "unknown"
    } else{
      shared.limit <- max(keep.scores) - shared_threshold
      final.scores <- keep.scores[keep.scores > shared.limit]
      final.names <- keep.names[keep.scores > shared.limit]
      if (length(final.scores)==1){
        classification <- final.names
        tissues <- final.names
      } else{
        classification <- "shared"
        tissues <- paste0(final.names,collapse = ",")
      }
    }
    build.df <- data.frame("SIGNAL"=row.df$SIGNAL,"classification"=classification,
                           "tissues"=tissues,stringsAsFactors = F)
    out.df <- rbind(out.df,build.df)
    setTxtProgressBar(pb,i)
  }
  return(out.df)
}
