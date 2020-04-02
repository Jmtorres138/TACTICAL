#' Calculate a tissue score vector for each SNP within a genetic credible set or for each index SNP
#'
#' This function takes the SNP annotations outputted from the annotate_snps function and the
#' provided annotation weigths to yield a set of tissue scores that indicate how much of the
#' fine-mapping/association signal can be apportioned to each tissue/cell type
#'
#' @param snp.annotated.df Dataframe of annotated SNP info, this is output from annotate_snps function
#' @param tissue_annotation_file Path to the input file of tissue annotation names and weights
#' @param genomic_annotation_file Path to the input file of genomic annotation names and weights]
#' @param ess.annot Optional: Name of genomic annotation requiring specificity scores (e.g. "coding")
#' @param ess.file Optional: Path to the specificity score file for ess.annot, (e.g. expression specificity scores for coding annotations)
#' @export
calculate_tissue_vectors <- function(snp.annotated.df,tissue_annotation_file,genomic_annotation_file,
                                     ess.annot=NULL,ess.file=NULL){
  "%&%" <- function(a,b) paste0(a,b) # just a shortcut for the paste function
  '%>%' <- magrittr::'%>%'
  tiss.annot.df <- data.table::fread(tissue_annotation_file)
  gen.annot.df <- data.table::fread(genomic_annotation_file)
  if (!is.null(ess.annot)){
    ess.df <- data.table::fread(ess.file)
  }
  out.df <- c()
  pb <- txtProgressBar(min=1,max=dim(snp.annotated.df)[1],style=3)
  for (r in 1:dim(snp.annotated.df)[1]){
    full.row.df <- snp.annotated.df[r,]
    row.df <- full.row.df %>% dplyr::select(.,-one_of("SIGNAL","SNPID","CHR","POS","VALUE"))
    # build annotation matrix for individual SNP
    name.vec <- names(row.df)
    tiss.vec <- purrr::map(name.vec,function(s){
      vec <- strsplit(x=s,split=".",fixed=T)[[1]]
      ifelse(length(vec)==2,vec[1],NA)
    }) %>% unlist(.) %>% na.omit(.) %>% as.character(.) %>% unique(.)  %>% sort(.)
    tiss.annot.vec <- purrr::map(name.vec,function(s){
      vec <- strsplit(x=s,split=".",fixed=T)[[1]]
      ifelse(length(vec)==2,vec[2],NA)
    }) %>% unlist(.) %>% na.omit(.) %>% as.character(.) %>% unique(.) %>% sort(.)
    gen.annot.vec <- purrr::map(name.vec,function(s){
      vec <- strsplit(x=s,split=".",fixed=T)[[1]]
      ifelse(length(vec)==1,vec[1],NA)
    }) %>% unlist(.) %>% na.omit(.) %>% as.character(.) %>% unique(.) %>% sort(.)
    annot.vec <- c(tiss.annot.vec,gen.annot.vec)
    annot.matrix <- matrix(nrow=length(annot.vec),ncol=length(tiss.vec))
    row.names(annot.matrix) <- annot.vec; colnames(annot.matrix) <- tiss.vec
    for (i in 1:length(annot.vec)){
      annot <- annot.vec[i]
      print.warning <- TRUE
      for (e in 1:length(tiss.vec)){
        tiss <- tiss.vec[e]
        aname <- ifelse(annot %in% tiss.annot.vec, tiss%&%"."%&%annot, annot)
        wgt <- ifelse(annot %in% tiss.annot.vec,
                      dplyr::filter(tiss.annot.df,V1==tiss,V2==annot)$V3,
                      dplyr::filter(gen.annot.df,V1==annot)$V2
        )
        if (!is.null(ess.annot)){
          if (aname==ess.annot){
            eval.snp <- dplyr::select(row.df,one_of(aname)) %>% as.numeric(.)
            if (eval.snp==1){
              snp.chr <- full.row.df$CHR; snp.pos <- full.row.df$POS
              sub.df <- dplyr::filter(ess.df,CHR==snp.chr,START<=snp.pos,END>=snp.pos)
              if (dim(sub.df)[1]==1){
                annot.matrix[i,e] <- (dplyr::select(sub.df,one_of(tiss)) %>% as.numeric(.)) * wgt
              } else if (dim(sub.df)[1]==0){
                if (print.warning==TRUE){
                  write("\nWARNING: SNP " %&% full.row.df$SNPID %&% " at SIGNAL " %&% full.row.df$SIGNAL %&%
                          " maps to annotation " %&% aname %&%
                          " but SNP does not map to feature in provided specificity file, please inspect",stdout())
                  print.warning<-FALSE
                }
                # when tissue specificity information is missing, defaults to equal value for each tissue
                default.val <- 1 / length(tiss.vec)
                annot.matrix[i,e] <- default.val * wgt
              } else if (dim(sub.df)[1]>1){
                if (print.warning==TRUE){
                  write("\nWARNING: SNP " %&% full.row.df$SNPID %&% " at SIGNAL " %&% full.row.df$SIGNAL %&%
                          " maps to annotation " %&% aname %&%
                          " but SNP maps to multiple features in provided specificity file, please inspect",stdout())
                  print.warning<-FALSE
                }
                # when tissue specificity information is ambiguous (i.e. snp maps to multiple features),
                # defaults to equal value for each tissue
                default.val <- 1 / length(tiss.vec)
                annot.matrix[i,e] <- default.val * wgt
              } else{
                write("\nError linking annotated SNP to specificity file, please inspect",stderr())
              }
            } else{
              annot.matrix[i,e] <- (dplyr::select(row.df,one_of(aname)) %>% as.numeric(.)) * wgt
            }
          } else{
            annot.matrix[i,e] <- (dplyr::select(row.df,one_of(aname)) %>% as.numeric(.)) * wgt
          }
        } else{
          annot.matrix[i,e] <- (dplyr::select(row.df,one_of(aname)) %>% as.numeric(.)) * wgt
        }

      }
    }
    # Sum and scale annotation vectors
    if (sum(colSums(annot.matrix)) > 0){
      tiss.wgts <- colSums(annot.matrix) / sum(colSums(annot.matrix))
    } else{
      tiss.wgts <- rep(0,length(tiss.vec))
      names(tiss.wgts) <- tiss.vec
    }
    # Obtain tissue scores for SNP by partitioning PPA value
    tiss.vector <- full.row.df$VALUE * tiss.wgts
    sub.df <- dplyr::select(full.row.df,one_of("SIGNAL","SNPID","CHR","POS","VALUE"))
    build.df <- cbind(sub.df,(as.data.frame(tiss.vector) %>% t(.)))
    out.df <- rbind(out.df,build.df)
    setTxtProgressBar(pb,r)
  }
  return(out.df)
}
