
  # spike in dinucleotide sequences
  # Load packages
  library(openxlsx)
  library(DESeq2)
  library(rtracklayer)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(BiocParallel)
  library(limma)
  library(apeglm)
  library(ashr)
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  
  
  # System specific parameters
  #setwd('/work/data/trad/')
  #register(MulticoreParam(12))
  setwd('~/Desktop/trad/')
  register(MulticoreParam(6))
  
  
  # Import lambda bed files
  bf <- dir('~/Desktop/tp53/bed/',pattern='spec2.bed.gz$',full.names=T)
  dat <- lapply(bf,import.bed)
  
  
  # Process
  tab <- lapply(dat,function(i) table(i$name))
  tab <- do.call(rbind,tab)
  row.names(tab) <- basename(bf)
  
  stderror <- function(x) sd(x)/sqrt(length(x))
  
  seqs <- colnames(tab)
  means <- apply(tab,2,mean)
  sem <- apply(tab,2,stderror)
  df <- data.frame(seqs=factor(seqs),means=means,sem=sem)
  
  ggplot(df,aes(x=seqs,y=means)) +
    geom_bar(stat='identity') +
    geom_errorbar(aes(ymin=means-sem,ymax=means+sem),width=0.2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  # Normalize by the sequence of the lambda genome
  gen <- readDNAStringSet('~/Desktop/tp53/lambda.fa.gz')
  bins <- slidingWindows(GRanges(seqnames='chrL',ranges=IRanges(start=1,width=48502)),width=2,step=1)
  
  
  lambdaseqs <- getSeq(gen,bins)
  tabls <- table(lambdaseqs[[1]])
  
  stopifnot(identical(names(means),names(tabls)))
  
  barplot(means / tabls,las=2)
  
  (means / tabls) / ((means / tabls)['CT'])
  
  
  
  
  