
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
  library(pheatmap)
  library(RColorBrewer)
  
  # System specific parameters
  #setwd('/work/data/trad/')
  #register(MulticoreParam(12))
  setwd('~/Desktop/trad/')
  register(MulticoreParam(6))
  
  
  
  # Import metadata
  coldata <- read.xlsx('bigwig_summary.xlsx')
  coldata[] <- lapply(coldata,factor)
  
  
  # Define parameters
  binSize <- 50e3 # sliding window width (in bp)
  binStep <- 50e3 # sliding window step (in bp)
  samples_to_plot <- c('IMR_100J_UVC','IMR_200J_UVB','Keratinocytes_200J_UVB','Melanocytes_200J_UVB','WI38_200J_UVB','TP53_200J_UVB','WI38RAF_200J_UVB')
  
  
  # Define genome bins
  human_chroms <- as(seqinfo(Hsapiens),'GRanges')[1:23]
  seqlevels(human_chroms) <- seqlevelsInUse(human_chroms)
  human_bins <- unlist(slidingWindows(human_chroms,width=binSize,step=binStep),use.names=F)
  names(human_bins) <- as.character(human_bins)
  lambda_chroms <- GRanges(seqnames='chrL',ranges=IRanges(start=1,end=48502))
  lambda_bins <- unlist(slidingWindows(lambda_chroms,width=20,step=10))
  names(lambda_bins) <- as.character(lambda_bins)
  
  
  # # Import bw files
  human_bw_files  <- BigWigFileList(dir('results/bw',full.names=T,pattern='_spec1.bw$'))
  lambda_bw_files <- BigWigFileList(dir('results/bw',full.names=T,pattern='_spec2.bw$'))
  #human_bw <- bplapply(human_bw_files,import,as='RleList')
  #human_bw <- lapply(human_bw,'[',seqlevels(human_bins))
  #lambda_bw <- bplapply(lambda_bw_files,import,as='RleList')
  #lambda_bw <- lapply(lambda_bw,'[',seqlevels(lambda_bins))
  #STOPPING POINT
  #saveRDS(human_bw,'human_bw.RDS')
  #saveRDS(lambda_bw,'lambda_bw.RDS')
  human_bw <- readRDS('human_bw.RDS')
  lambda_bw <- readRDS('lambda_bw.RDS')
  
  
  # Compute counts
  lambda_counts <- bplapply(lambda_bw,function(i) {
    as.integer(binnedAverage(bins=lambda_bins,numvar=i,varname='x')$x * width(lambda_bins))
  })
  lambda_counts <- do.call(cbind,lambda_counts)
  colnames(lambda_counts) <- unlist(strsplit(basename(lambda_bw_files),'_possort_markdup_spec2.bw'))
  row.names(lambda_counts) <- names(lambda_bins)
  human_counts <- bplapply(human_bw,function(i) {
    as.integer(binnedAverage(bins=human_bins,numvar=i,varname='x')$x * width(human_bins))
  })
  human_counts <- do.call(cbind,human_counts)
  colnames(human_counts) <- unlist(strsplit(basename(human_bw_files),'_possort_markdup_spec1.bw'))
  row.names(human_counts) <- names(human_bins)
  
  
  # Construct SummarizedExperiment object
  # sanity check
  stopifnot(colnames(human_counts) == coldata$library_ID)
  stopifnot(colnames(lambda_counts) == coldata$library_ID)
  counts <- rbind(human_counts,lambda_counts)
  bins <- suppressWarnings(c(human_bins,lambda_bins))
  se <- SummarizedExperiment(assays=list(counts=counts),
                             rowRanges=bins,
                             colData=coldata)
  # STOPPING POINT
  # # saveRDS(se,'se_020623.RDS')
  # # se <- readRDS('se_020623.RDS')
  
  
  # Estimate Size Factors from spike-in
  dds <- DESeqDataSet(se,design=~batch+condition)
  dds$condition <- relevel(dds$condition,ref='IMR_200J_UVB') # set reference level
  controls <- grepl('^chrL:',row.names(se))
  dds <- estimateSizeFactors(dds,type='ratio',locfunc=median,controlGenes=controls) # can be ratio, poscounts, or iterate
  
  
  # Differential susceptibility analysis
  dds <- dds[!grepl('^chrL:',row.names(se))] # remove spike-in counts
  dds <- DESeq(dds,parallel=TRUE)
  Contrast <- c('condition','TP53_200J_UVB','IMR_200J_UVB')
  #Contrast <- c('condition','IMR_100J_UVC','IMR_200J_UVB')
  #Contrast <- c('condition','WI38_200J_UVB','IMR_200J_UVB')
  #Contrast <- c('condition','KMT5BC_200J_UVB','IMR_200J_UVB')
  res <- results(dds,contrast=Contrast)
  sum(res$padj < 0.05,na.rm=TRUE)
  sum(res$padj < 0.05,na.rm=TRUE) / nrow(res)
  sum(res$padj < 0.05 & res$log2FoldChange>1,na.rm=T)
  sum(res$padj < 0.05 & res$log2FoldChange<0,na.rm=T)
  sum(res$padj < 0.00001,na.rm=TRUE)
  
  resASH <- lfcShrink(dds, contrast=Contrast, type='ashr') # log fold-change shrinkage for visualizations
  #resASH$lambda_control <- grepl('^chrL:',row.names(resASH))
  #resASH <- resASH[complete.cases(resASH),] # filtering
  #DSR_thresh <- c(0.05,1)
  #resASH$DSR <- ifelse(resASH$padj < DSR_thresh[1] & abs(resASH$log2FoldChange) > DSR_thresh[2],TRUE,FALSE)
  #resASH$minus_log10_pval <- -1*log10(resASH$pvalue)
  
  #saveRDS(resASH,file='TP53_resASH.RDS')
  
  
  # PCA
  dds <- dds[,dds$condition %in% samples_to_plot]
  dds$condition <- droplevels(dds$condition)
  rld <- rlog(dds, blind=TRUE)
  mat <- assay(rld)
  mm <- model.matrix(~condition,colData(rld))
  mat <- limma::removeBatchEffect(mat,batch=rld$batch,design=mm)
  assay(rld) <- mat
  ntop <- round(0.1*nrow(dds)) # top 10%
  plotPCA(rld,ntop=ntop)
  
  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)