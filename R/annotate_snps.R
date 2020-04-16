#' Annotate SNPs from genetic fine-mapping or GWAS index variants
#'
#' This function annotates SNPs to genomic annotations (e.g. coding sequence) and
#' tissue-level annotations (e.g. chromatin states). Each SNP is mapped against each
#' set of annotations contained within BED files listed within the genomic annotation file and
#' tissue annotation file. The annotation names listed in the annotation files must be consistent
#' with the annotation names listed within the BED files.
#'
#' @param snp_file Path to the input SNP file
#' @param tissue_path_file Path to the input file of paths to tissue annotation bed files
#' @param tissue_annotation_file Path to the input file of tissue annotation names and weights
#' @param genomic_path_file Path to the input file of paths to genomic annotation bed files
#' @param genomic_annotation_file Path to the input file of genomic annotation names and weights
#' @export
annotate_snps <- function(snp_file,tissue_path_file,tissue_annotation_file,
                          genomic_path_file,genomic_annotation_file){
  "%&%" <- function(a,b) paste0(a,b) # just a shortcut for the paste function
  '%>%' <- magrittr::'%>%'
  '%over%' <- IRanges::'%over%'
  snp.df <- data.table::fread(snp_file,header=T)
  tiss.df <- data.table::fread(file=tissue_path_file,header=F)
  tiss.annot.df <- data.table::fread(file=tissue_annotation_file,header=F)
  gen.df <- data.table::fread(file=genomic_path_file,header=F)
  gen.annot.df <- data.table::fread(file=genomic_annotation_file)
  out.df <- snp.df
  snp.gr <- GenomicRanges::GRanges(seqnames=snp.df$CHR,IRanges::IRanges(start=snp.df$POS,end=snp.df$POS))
  len <- dim(tiss.df)[1] * length(unique(tiss.annot.df$V2))
  pb <- txtProgressBar(min=0,max=len,style=3)
  count <- 0
  for (i in 1:dim(tiss.df)[1]){
    tissue <- tiss.df$V1[i]; pth <- tiss.df$V2[i]
    bed.df <- data.table::fread(file=pth,header=F,stringsAsFactors = F,sep="\t")
    for (a in unique(tiss.annot.df$V2)){
      count <- count + 1
      aname <- tissue%&%"."%&%a
      sub.df <- dplyr::filter(bed.df,V4==a)
      annot.gr <- GenomicRanges::GRanges(seqnames=sub.df$V1,IRanges::IRanges(start=sub.df$V2+1,end=sub.df$V3+1)) 
      eval.vec <- (snp.gr %over% annot.gr) %>% as.integer(.)
      out.df <- cbind(out.df,eval.vec)
      names(out.df)[dim(out.df)[2]] <- aname
      setTxtProgressBar(pb,count)
    }
  }
  for (i in 1:dim(gen.df)[1]){
    gannot <- gen.df$V1[i]; pth <- gen.df$V2[i]
    bed.df <- data.table::fread(file=pth,header=F,stringsAsFactors = F,,sep="\t")
    for (a in unique(gen.annot.df$V1)){
      count <- count + 1
      sub.df <- dplyr::filter(bed.df,V4==a)
      annot.gr <- GenomicRanges::GRanges(seqnames=sub.df$V1,IRanges::IRanges(start=sub.df$V2+1,end=sub.df$V3+1))
      eval.vec <- (snp.gr %over% annot.gr) %>% as.integer(.)
      out.df <- cbind(out.df,eval.vec)
      names(out.df)[dim(out.df)[2]] <- gannot
      setTxtProgressBar(pb,count)
    }
  }
  return(out.df)
}
