# TACTICAL: Tissue of ACTion scores for Investigating Complex trait Associations at Loci 

This package takes genetic associations or fine-mapped genetic credible sets and systemically integrates them with functional annotations to obtain *tissue of action* (TOA) scores that can guide mechanistic investigation of genetic signals from genome-wide association studies.  

### Prerequisites 

The package was developed in R (version 3.6.0) and requires the following R packages: 

* data.table (1.12.8)
* dplyr (0.8.4)
* GenomicRanges (1.36.1)
* devtools (2.2.2)
    * __`install_github`__ function used for package installation
* **Note**: functions from ggplot2 (3.3.0) and gridExtra (2.3) are used to plot PCA output in this tutorial 
    

### Installation 

The package can be directly installed from [GitHub](https://github.com/Jmtorres138/TACTICAL) using this R command: 

```
devtools::install_github("jmtorres138/TACTICAL")
```

### Suggested usage 

TACTICAL provides a simple way to systematically integrate genetic information from trait-associated loci with functional genomic annotations and gene expression to obtain *tissue-of-action* (TOA) scores. These scores can be used to:

* classify signals to most likely tissues of action. 
* prioritise genetic signals for experimental validation in a particular cell or tissue type.
* guide the identification of causal gene(s) - at specific loci - by informing which tissue(s) are most appropriate for analyses such as eQTL colocalisation or integration with chromatin conformation capture (3C) approaches.
* inform *process-specific* polygenic risk scores. 

TACTICAL is simple to run and flexible; you need only provide the input data you are interested in evaluating. The data accepted include: 

* `Genetic information`: Credible sets from bayesian fine-mapping or index variants from GWAS. If using index SNPs, we recommended using conditionally independent SNPs or LD-pruned SNPs. 
* `Genomic annotations`: Any type of interval data that can be formatted in a BED file. This can include chromatin segmentation states, ChIP-seq peaks, chromatin accessible regions (i.e. DHS or ATAC-seq peaks), lowly/highly methylated regions, coding sequence, untranslated regions, etc.
* `Expression specificity scores`: When comparing across a set of tissues or cell types, you can provide expression specificity scores (described below) that inform the extent of tissue-specific gene expression in each evaluated tissue. 

Although there are many possible data combinations the user may wish to explore, we offer a few suggestions: 

1. As some annotations are more enriched for genome-wide significant SNPs (or SNP heritability), you may consider explicitly using these enrichment values as weights within TACTICAL. Although this package does not estimate enrichment *per se*, there are many software programs available that do, such as [fgwas](https://github.com/joepickrell/fgwas), [GARFIELD](https://www.bioconductor.org/packages/release/bioc/html/garfield.html), and [GoShifter](https://github.com/immunogenomics/goshifter).
2. If you are interested in prioritising genetic signals within a particular tissue, then we suggest using as many relevant annotations as available (i.e. from a variety of molecular assays); though you may not want to include "repressed" or low-signal annotations in your tissue or cell-type of interest. 
3. On the other hand, if you are interested in profiling signals across a range of candidate tissues, then we strongly recommend restricting the provided annotations to those are available for all evaluated tissues. For example, you may not want to use ATAC-seq peaks if they are only available for a subset of tissues.  
