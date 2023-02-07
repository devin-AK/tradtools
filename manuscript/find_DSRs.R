  
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
  
  
  # System specific parameters
  #setwd('/work/data/trad/')
  #register(MulticoreParam(12))
  setwd('~/Desktop/trad/')
  register(MulticoreParam(6))
  
  
  # Import metadata
  coldata <- read.xlsx('bigwig_summary.xlsx')
  coldata[] <- lapply(coldata,factor)

  
  # Define genome bins
  human_chroms <- as(seqinfo(Hsapiens),'GRanges')[1:23]
  seqlevels(human_chroms) <- seqlevelsInUse(human_chroms)
  human_bins <- unlist(slidingWindows(human_chroms,width=100e3,step=50e3),use.names=F)
  names(human_bins) <- as.character(human_bins)
  lambda_chroms <- GRanges(seqnames='chrL',ranges=IRanges(start=1,end=48502))
  lambda_bins <- unlist(slidingWindows(lambda_chroms,width=20,step=10))
  names(lambda_bins) <- as.character(lambda_bins)
  
  
  # Import bw files
  human_bw_files  <- BigWigFileList(dir('results/bw',full.names=T,pattern='_spec1.bw$'))
  human_bw <- bplapply(human_bw_files,import,as='RleList')
  human_bw <- lapply(human_bw,'[',seqlevels(human_bins))
  lambda_bw_files <- BigWigFileList(dir('results/bw',full.names=T,pattern='_spec2.bw$'))
  lambda_bw <- bplapply(lambda_bw_files,import,as='RleList')
  lambda_bw <- lapply(lambda_bw,'[',seqlevels(lambda_bins))
  # STOPPING POINT
  # saveRDS(human_bw,'human_bw.RDS')
  # saveRDS(lambda_bw,'lambda_bw.RDS')
  # human_bw <- readRDS('human_bw.RDS')
  # lambda_bw <- readRDS('lambda_bw.RDS')
  
  
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
  
  
  # Estimate Size Factors (to account for differences in sequencing depth)
  dds <- DESeqDataSet(se,design=~batch+condition)
  dds$condition <- relevel(dds$condition,ref='IMR_200J_UVB') # set reference level
  controls <- grepl('^chrL:',row.names(se))
  dds <- estimateSizeFactors(dds,type='ratio',locfunc=median,controlGenes=controls) # can be ratio, poscounts, or iterate
  
  
  # Differential susceptibility analysis
  dds <- DESeq(dds,parallel=TRUE)
  
  
  # Results
  #Contrast <- c('condition','IMR_100J_UVC','IMR_200J_UVB')
  Contrast <- c('condition','TP53_200J_UVB','IMR_200J_UVB')
  res <- results(dds,contrast=Contrast)
  resASH <- lfcShrink(dds, contrast=Contrast, type='ashr') # log fold-change shrinkage for visualizations
  resASH$lambda_control <- controls
  #resASH <- resASH[complete.cases(resASH),] # filtering
  DSR_thresh <- c(0.00001,1)
  resASH$DSR <- ifelse(resASH$padj < DSR_thresh[1] & abs(resASH$log2FoldChange) > DSR_thresh[2],TRUE,FALSE)
  resASH$minus_log10_pval <- -1*log10(resASH$pvalue)
  
  
  # Diagnostic Plots
  DESeq2::plotDispEsts(dds)
  DESeq2::plotSparsity(dds)
  DESeq2::plotMA(res,ylim=c(-4,4))
  DESeq2::plotMA(resASH,ylim=c(-4,4))
  
  
  # Spike in plots
  ps1 <- ggplot(as.data.frame(resASH),aes(x=baseMean,y=log2FoldChange)) +
    geom_point(aes(color=lambda_control)) +
    geom_hline(yintercept=c(-1,1),linetype='dashed') +
    scale_color_manual(values=c('gray','#9ECAE1')) +
    scale_x_continuous(trans='log10',limits=c(0.01,15000)) +
    ylim(c(-4,4)) +
    labs(color='lambda control') +
    xlab('mean susceptibility') +
    ylab('log fold-change') +
    guides(color='none') +
    theme_bw()
  ps2 <- ggplot(as.data.frame(resASH),aes(x=log2FoldChange,y=minus_log10_pval)) +
    geom_point(aes(color=lambda_control)) +
    scale_color_manual(values=c('gray','#9ECAE1')) +
    xlim(c(-4,4)) +
    ylim(c(0,20)) +
    ylab('-log10(pvalue)') +
    geom_hline(yintercept = -log10(DSR_thresh[1]), linetype='dashed') +
    geom_vline(xintercept = c(-1*DSR_thresh[2],DSR_thresh[2]), linetype='dashed') + 
    theme_bw()
  spike_plots <- ggarrange(ps1,ps2,ncol=2,widths=c(1,1.5))
  annotate_figure(spike_plots,top=paste(Contrast[2],'vs',Contrast[3]))
  
  
  # MA and Volcano plots
  pm1 <- ggplot(as.data.frame(resASH),aes(x=baseMean,y=log2FoldChange)) +
    geom_point(aes(color=DSR)) +
    geom_hline(yintercept=c(-1,1),linetype='dashed') +
    scale_color_manual(values=c('gray','#FC9272')) +
    scale_x_continuous(trans='log10',limits=c(0.01,15000)) +
    ylim(c(-4,4)) +
    labs(color='lambda control') +
    xlab('mean susceptibility') +
    ylab('log fold-change') +
    guides(color='none') +
    theme_bw()
  pv2 <- ggplot(as.data.frame(resASH),aes(x=log2FoldChange,y=minus_log10_pval)) +
    geom_point(aes(color=DSR)) +
    #geom_point(aes(fill=lambda_control,color=DSR),pch=21) +
    scale_color_manual(values=c('gray','#FC9272')) +
    #scale_fill_manual(values=c('gray','dodgerblue')) +
    xlim(c(-4,4)) +
    ylim(c(0,20)) +
    ylab('-log10(pvalue)') +
    geom_hline(yintercept = -log10(DSR_thresh[1]), linetype='dashed') +
    geom_vline(xintercept = c(-1*DSR_thresh[2],DSR_thresh[2]), linetype='dashed') + 
    theme_bw()
  DSR_plots <- ggarrange(pm1,pv2,ncol=2,widths=c(1,1.5))
  annotate_figure(DSR_plots,top=paste(Contrast[2],'vs',Contrast[3]))
  
  
  # MA and Volcano plots option 2
  pm1 <- ggplot(as.data.frame(resASH),aes(x=baseMean,y=log2FoldChange)) +
    geom_point(aes(fill=lambda_control,color=DSR),pch=21) +
    scale_color_manual(values=c('gray','#FC9272')) +
    scale_fill_manual(values=c('gray','dodgerblue')) +
    geom_hline(yintercept=c(-1,1),linetype='dashed') +
    scale_x_continuous(trans='log10',limits=c(0.01,15000)) +
    ylim(c(-4,4)) +
    labs(color='lambda control') +
    xlab('mean susceptibility') +
    ylab('log fold-change') +
    guides(color='none') +
    theme_bw()
  pv2 <- ggplot(as.data.frame(resASH),aes(x=log2FoldChange,y=minus_log10_pval)) +
    #geom_point(aes(color=DSR)) +
    geom_point(aes(fill=lambda_control,color=DSR),pch=21) +
    scale_color_manual(values=c('gray','#FC9272')) +
    scale_fill_manual(values=c('gray','dodgerblue')) +
    xlim(c(-4,4)) +
    ylim(c(0,20)) +
    ylab('-log10(pvalue)') +
    geom_hline(yintercept = -log10(DSR_thresh[1]), linetype='dashed') +
    geom_vline(xintercept = c(-1*DSR_thresh[2],DSR_thresh[2]), linetype='dashed') + 
    theme_bw()
  DSR_plots <- ggarrange(pm1,pv2,ncol=2,widths=c(1,1.5))
  annotate_figure(DSR_plots,top=paste(Contrast[2],'vs',Contrast[3]))
  
  
  # PCA plots
  vsd <- vst(dds,blind=T,fitType='mean')
  p1 <- plotPCA(vsd,'condition')
  p2 <- plotPCA(vsd,'batch')
  mat <- assay(vsd)
  mat <- limma::removeBatchEffect(mat,vsd$batch)
  assay(vsd) <- mat
  p3 <- plotPCA(vsd,'condition')
  p4 <- plotPCA(vsd,'batch')
  grid.arrange(p1,p2,p3,p4,nrow=2)
  
  
  # DSR calculation
  dds <- DESeq(dds,fitType='mean')
  res <- results(dds,contrast=c('condition','TP53_200J_UVB','IMR_200J_UVB'))
  
  
  EnhancedVolcano(res,lab=rownames(res),x='log2FoldChange',y='pvalue')
  
  
  # subset (optional)
  # samples_to_keep <- c('IMR_100J_UVC','IMR_200J_UVB','Keratinocytes_200J_UVB','Melanocytes_200J_UVB','WI38_200J_UVB','TP53_200J_UVB')
  # dds <- dds[,dds$condition %in% samples_to_keep]
  # dds <- estimateSizeFactors(dds,type='iterate',controlGenes=controlGenes)