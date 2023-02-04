
  # Load packages
  library(openxlsx)
  library(DESeq2)
  library(rtracklayer)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(BiocParallel)
  library(limma)

  
  # System specific parameters
  setwd('/work/data/trad/')
  register(MulticoreParam(12))
  
  
  # Import metadata
  metadata <- read.xlsx('bigwig_summary.xlsx')
  
  
  # Define genome bins
  human_chroms <- as(seqinfo(Hsapiens),'GRanges')[1:23]
  seqlevels(human_chroms) <- seqlevelsInUse(human_chroms)
  human_bins <- unlist(slidingWindows(human_chroms,width=100e3,step=50e3),use.names=F)
  names(human_bins) <- as.character(human_bins)
  lambda_chroms <- GRanges(seqnames='chrL',ranges=IRanges(start=1,end=48502))
  lambda_bins <- unlist(slidingWindows(lambda_chroms,width=1000,step=500))
  names(lambda_bins) <- as.character(lambda_bins)
  
  
  # Import bw files
  human_bw_files  <- BigWigFileList(dir('results/bw',full.names=T,pattern='_spec1.bw$'))
  human_bw <- bplapply(human_bw_files,import,as='RleList')
  human_bw <- lapply(human_bw,'[',seqlevels(human_bins))
  lambda_bw_files <- BigWigFileList(dir('results/bw',full.names=T,pattern='_spec2.bw$'))
  lambda_bw <- bplapply(lambda_bw_files,import,as='RleList')
  lambda_bw <- lapply(lambda_bw,'[',seqlevels(lambda_bins))
  
  
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
  
  
  # Construct SummarizedExperiment and DESeqDataSet objects
  # sanity check
  identical(colnames(human_counts),metadata$library_ID)
  identical(colnames(lambda_counts),metadata$library_ID)
  counts <- rbind(human_counts,lambda_counts)
  bins <- c(human_bins,lambda_bins)
  se <- SummarizedExperiment(assays=list(counts=counts),
                             rowRanges=bins,
                             colData=metadata)
  
  # 
  # hc <- SummarizedExperiment(assays=list(counts=human_counts),
  #                            rowRanges=human_bins,
  #                            colData=metadata)
  # lc <- SummarizedExperiment(assays=list(counts=lambda_counts),
  #                            rowRanges=lambda_bins,
  #                            colData=metadata)
  # hdds <- DESeqDataSet(hc,design=~batch+sample)
  # ldds <- DESeqDataSet(lc,design=~batch+sample)
  
  
  # Estimate size factors from spike in
  controlGenes <- grep('^chrL:',row.names(counts),value=T)
  dds <- DESeqDataSet(se,design=~batch+condition)
  dds$condition <- relevel(dds$condition,ref='IMR_200J_UVB')
  dds <- estimateSizeFactors(dds,type='iterate',controlGenes=controlGenes)
  # dds <- readRDS('dds_020323.RDS')
  
  
  # subset (optional)
  samples_to_keep <- c('IMR_100J_UVC','IMR_200J_UVB','Keratinocytes_200J_UVB','Melanocytes_200J_UVB','WI38_200J_UVB','TP53_200J_UVB')
  dds <- dds[,dds$condition %in% samples_to_keep]
  dds <- estimateSizeFactors(dds,type='iterate',controlGenes=controlGenes)
  
  
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
  
  
  
  tmp <- res[!grepl('^chrL:',row.names(res)),]
  
  EnhancedVolcano(tmp,lab=rownames(tmp),x='log2FoldChange',y='pvalue')
  
  
  
  dds <- readRDS('dds_020323.RDS')
  samples_to_keep <- c('IMR_100J_UVC','IMR_200J_UVB','Keratinocytes_200J_UVB','Melanocytes_200J_UVB','WI38_200J_UVB','TP53_200J_UVB')
  samples_to_keep <- c(samples_to_keep,'WI38RAF_200J_UVB')
  dds <- dds[,dds$condition %in% samples_to_keep]
  design(dds) <- formula(~batch + condition)
  dds$condition <- droplevels(dds$condition)
  lc <- dds[grepl('^chrL:',row.names(dds)),]
  lc <- estimateSizeFactors(lc)
  sizeFactors(dds) <- sizeFactors(lc)
  dds <- DESeq(dds)
  res <- results(dds,contrast=c('condition','TP53_200J_UVB','IMR_200J_UVB'))
  EnhancedVolcano(res[!grepl('^chrL:',row.names(res)),],title='TP53',lab=NA,x='log2FoldChange',y='pvalue')
  EnhancedVolcano(res[grepl('^chrL:',row.names(res)),],title='lambda',lab=NA,x='log2FoldChange',y='pvalue')
  
  
  EnhancedVolcano(res[!grepl('^chrL:',row.names(res)),],x='log2FoldChange',y='pvalue',lab=NA)
  
  spike <- grep('^chrL:',row.names(res),value=T)
  EnhancedVolcano(res,lab=row.names(res),x='log2FoldChange',y='padj')
  
  EnhancedVolcano(res,lab=rownames(res),x='log2FoldChange',y='pvalue',colCustom=spike)
  
  