
  # Define parameters
  rm(list=ls())
  gc()
  samples_to_plot <- c('WI38_200J_UVB','IMR_200J_UVB','Keratinocytes_200J_UVB','Melanocytes_200J_UVB','IMR_100J_UVC','TP53_200J_UVB')
  colors <- setNames(brewer.pal(n=length(samples_to_plot),'Paired'),samples_to_plot)
  #Contrast <- c('condition','WI38_200J_UVB','IMR_200J_UVB')
  #Contrast <- c('condition','Keratinocytes_200J_UVB','IMR_200J_UVB')
  #Contrast <- c('condition','Melanocytes_200J_UVB','IMR_200J_UVB')
  #Contrast <- c('condition','IMR_100J_UVC','IMR_200J_UVB')
  Contrast <- c('condition','TP53_200J_UVB','IMR_200J_UVB')
  pre_filter_threshold <- 100
  DSR_thresh <- c(0.00001,1)
  include_spike_reads <- FALSE
  binSize <- 50e3 # sliding window width (in bp)
  binStep <- 50e3 # sliding window step (in bp)



  # Load packages
  library(openxlsx)
  library(DESeq2)
  library(rtracklayer)
  library(GenomicRanges)
  library(circlize)
  library(data.table)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(BiocParallel)
  library(limma)
  library(apeglm)
  library(ashr)
  library(ggplot2)
  library(ggpubr)
  library(pheatmap)
  library(RColorBrewer)
  library(factoextra)
  library(NbClust)
  library(ggrastr)
  
  
  # Define functions
  count_DSRs_in_tiles <- function(dsr, tiles, padj_thresh=0.05, log2fc_thresh=1) {
    stopifnot(is(dsr,'GRanges'))
    stopifnot(is(tiles,'GRanges'))
    stopifnot(!is.null(dsr$padj))
    stopifnot(!is.null(dsr$log2FoldChange))
    # process
    dsr.up <- dsr[dsr$padj < padj_thresh & dsr$log2FoldChange > abs(log2fc_thresh)]
    dsr.dn <- dsr[dsr$padj < padj_thresh & dsr$log2FoldChange < -1*(abs(log2fc_thresh))]
    counts.up <- countOverlaps(tiles,dsr.up,type='any',ignore.strand=TRUE)
    counts.dn <- countOverlaps(tiles,dsr.dn,type='any',ignore.strand=TRUE)
    bed <- data.frame(chr=seqnames(tiles),start=start(tiles),end=end(tiles),up=counts.up,dn=counts.dn)
    return(bed)
  }
  
  
  
  DSR_plot3 <- function(bed, col_range=c(0,30), ylim=col_range, outline_points=F, hline_interval=5, cex_modifier=1, bg.border=NA) {
    require(circlize)
    require(RColorBrewer)
    # Initialize
    circos.par(cell.padding = c(0, 0, 0, 0),gap.after=2,start.degree=90,'track.height'=0.3)
    circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22),plotType=c('ideogram','labels'),axis.labels.cex=cex_modifier*0.4,labels.cex=cex_modifier*0.8)
    # up
    bed_up <- bed[,1:4]
    color_ramp_up <- brewer.pal(9,'Reds')[3:7]
    col_breaks_up <- seq(col_range[1],col_range[2],length.out=length(color_ramp_up)+1)
    col_breaks_up[1] <- min(min(bed_up[,4]),col_range[1])
    col_breaks_up[length(col_breaks_up)] <- max(max(bed_up[,4]),col_range[2])
    circos.genomicTrack(bed_up, panel.fun = function(region, value, ...) {
      for(h in seq(min(ylim),max(ylim),by=hline_interval)) {
        circos.lines(CELL_META$cell.xlim,c(h,h),lty=2,lwd=0.25,col='#AAAAAA')
      }
      v <- value[[1]]
      minv <- min(v)
      maxv <- max(v)
      cex <- (v - minv)/(maxv - minv)
      cex <- cex * cex_modifier
      breaks_up <- cut(v,breaks=col_breaks_up,include.lowest=T)
      names(color_ramp_up) <- levels(breaks_up)
      cols_up <- unname(color_ramp_up[breaks_up])
      if(outline_points) circos.genomicPoints(region, value, cex=cex, pch = 21, bg = cols_up, ...) else {
        circos.genomicPoints(region, value, cex=cex, pch = 16, col = cols_up, ...)
      }
    }, ylim = ylim, bg.border=bg.border, bg.col=add_transparency('#DE2D26',0.95))
    # middle border
    circos.track(ylim = c(0, 1), track.height = 0.01, bg.border=NA, bg.col='gray40')
    # dn
    bed_dn <- cbind(bed[,1:3],bed[,5])
    color_ramp_dn <- brewer.pal(9,'Blues')[3:7]
    col_breaks_dn <- seq(col_range[1],col_range[2],length.out=length(color_ramp_dn)+1)
    col_breaks_dn[1] <- min(min(bed_dn[,4]),col_range[1])
    col_breaks_dn[length(col_breaks_dn)] <- max(max(bed_dn[,4]),col_range[2])
    circos.genomicTrack(bed_dn, panel.fun = function(region, value, ...) {
      for(h in -1*seq(min(ylim),max(ylim),by=hline_interval)) {
        circos.lines(CELL_META$cell.xlim,c(h,h),lty=2,lwd=0.25,col='#AAAAAA')
      }
      v <- value[[1]]
      minv <- min(v)
      maxv <- max(v)
      cex <- (v - minv)/(maxv - minv)
      cex <- cex * cex_modifier
      breaks_dn <- cut(v,breaks=col_breaks_dn,include.lowest=T)
      names(color_ramp_dn) <- levels(breaks_dn)
      cols_dn <- unname(color_ramp_dn[breaks_dn])
      v_dn <- v*-1
      if(outline_points) circos.genomicPoints(region, v_dn, cex=cex, pch = 21, bg = cols_dn, ...) else {
        circos.genomicPoints(region, v_dn, cex=cex, pch = 16, col = cols_dn, ...)
      }
    }, ylim = rev(ylim)*-1, bg.border=bg.border, bg.col=add_transparency('#3182BD',0.95))
    circos.clear()
  }
  
  
  
  
  plotDispEsts2 <- function(dds){
    require(dplyr)
    require(reshape2)
    require(ggplot2)
    as.data.frame(mcols(dds)) %>% 
      select(baseMean, dispGeneEst, dispFit, dispersion) %>% 
      melt(id.vars="baseMean") %>% 
      filter(baseMean>0) %>% 
      ggplot(aes(x=baseMean, y=value, colour=variable)) + 
      rasterize(geom_point(size=0.1),dpi=600) +
      scale_x_log10() + 
      scale_y_log10() + 
      theme_bw() + 
      ylab("Dispersion") + 
      xlab("BaseMean") +
      scale_colour_manual(
        values=c("Black", "#e41a1c", "#377eb8"), 
        breaks=c("dispGeneEst", "dispFit", "dispersion"), 
        labels=c("Estimate", "Fit", "Final"),
        name=""
      ) +
      guides(colour = guide_legend(override.aes = list(size=2)))
  }
  
  
  
  remove_geoms <- function(x, geom_type) {
    # Find layers that match the requested type.
    selector <- sapply(x$layers,
                       function(y) {
                         class(y$geom)[1] == geom_type
                       })
    # Delete the layers.
    x$layers[selector] <- NULL
    x
  }
  
  
  
  
  # System specific parameters
  #setwd('/work/data/trad/')
  #register(MulticoreParam(12))
  setwd('~/Desktop/trad/')
  register(MulticoreParam(6))
  
  
  
  # Import metadata
  coldata <- read.xlsx('bigwig_summary.xlsx')
  coldata[] <- lapply(coldata,factor)
  coldata$condition <- factor(coldata$condition,levels=samples_to_plot)
  coldata$condition <- droplevels(coldata$condition)

  
  
  # Define genome bins
  human_chroms <- as(seqinfo(Hsapiens),'GRanges')[1:23]
  seqlevels(human_chroms) <- seqlevelsInUse(human_chroms)
  human_bins <- unlist(slidingWindows(human_chroms,width=binSize,step=binStep),use.names=F)
  names(human_bins) <- as.character(human_bins)
  lambda_chroms <- GRanges(seqnames='chrL',ranges=IRanges(start=1,end=48502))
  lambda_bins <- unlist(slidingWindows(lambda_chroms,width=20,step=10))
  names(lambda_bins) <- as.character(lambda_bins)
  
  
  
  # Import bw files
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
  se <- se[,se$condition %in% samples_to_plot]
  # STOPPING POINT
  # # saveRDS(se,'se_020623.RDS')
  # # se <- readRDS('se_020623.RDS')
  
  
  
  # Estimate Size Factors from spike-in
  dds <- DESeqDataSet(se,design=~batch+condition)
  dds$condition <- relevel(dds$condition,ref='IMR_200J_UVB') # set reference level
  controls <- grepl('^chrL:',row.names(se))
  dds <- estimateSizeFactors(dds,type='ratio',locfunc=median,controlGenes=controls)
  
  
  
  # Differential susceptibility analysis
  if(!include_spike_reads) dds <- dds[!grepl('^chrL:',row.names(se))] # remove spike-in counts
  dds <- dds[rowSums(counts(dds)) >= pre_filter_threshold] # pre-filter
  dds <- DESeq(dds,parallel=TRUE)
  dds <- dds[,dds$condition %in% samples_to_plot]
  #dds$condition <- droplevels(dds$condition)
  
  
  
  # Dispersion estimates
  # -------------------------------------------------------------------------- #
  #pdf(file='figures/DSRs/disp_est.pdf',width=4,height=2.4)
  plotDispEsts2(dds)
  #dev.off()
  # -------------------------------------------------------------------------- #
  
  
  
  # PCA
  rld <- rlog(dds, blind=TRUE)
  mat <- assay(rld)
  mm <- model.matrix(~condition,colData(rld))
  mat <- limma::removeBatchEffect(mat,batch=rld$batch,design=mm)
  assay(rld) <- mat
  ntop <- round(0.1*nrow(dds)) # top 10%
  p.pca <- plotPCA(rld,ntop=ntop) +
          ggtitle(paste0('PCA (bin size: ',binSize)) +
          scale_color_manual(values=colors) +
          theme_classic(base_size=6)
  p.pca <- remove_geoms(p.pca,'GeomPoint') + geom_point(size=1)
  p.pca
  LEGEND <- get_legend(p.pca)
  # -------------------------------------------------------------------------- #
  #pdf(file='figures/DSRs/LEGEND.pdf',width=4,height=4)
  as_ggplot(LEGEND)
  #dev.off()
  # -------------------------------------------------------------------------- #
  
  
  
  
  # Cluster analysis
  rv <- rowVars(assay(rld))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
  pca <- prcomp(t(assay(rld)[select, ]))
  #p.clus <- fviz_eig(pca) +
  #  theme_minimal(base_size=7)
  p.clus <- fviz_eig(pca,geom='line',size=0.2,linetype='dashed') +
    theme_classic(base_size=6) +
    scale_y_log10()
  p.clus <- remove_geoms(p.clus,'GeomPoint')
  #saveRDS(layer_data(p.clus),file=paste0('figures/DSRs/pclus_binSize_',binSize,'.RDS'))
  #saveRDS(se,file=paste0('figures/DSRs/se_binSize_',binSize,'.RDS'))
  # -------------------------------------------------------------------------- #
  #pdf(file=paste0('figures/DSRs/binSize_',binSize,'.pdf'),width=3,height=0.8)
    ggarrange(p.pca + theme(legend.position='none'),p.clus,ncol=2,widths=c(2,1))
    #ggarrange(p.pca + theme(legend.position='none'),p.clus,ncol=2,widths=c(2,1),align='hv')
  #dev.off()
  # -------------------------------------------------------------------------- #
  
  
  
  # DSR results
  res <- results(dds,contrast=Contrast)
  resASH <- lfcShrink(dds, contrast=Contrast, type='ashr') # log fold-change shrinkage for visualizations
  resASH$lambda_control <- grepl('^chrL:',row.names(resASH))
  #resASH <- resASH[complete.cases(resASH),] # filtering
  resASH$DSR <- ifelse(resASH$padj < DSR_thresh[1] & abs(resASH$log2FoldChange) > DSR_thresh[2],TRUE,FALSE)
  resASH$minus_log10_pval <- -1*log10(resASH$pvalue)
  
  
  
  # MA and Volcano plots
  pm1 <- ggplot(as.data.frame(resASH),aes(x=baseMean,y=log2FoldChange)) +
    rasterize(geom_point(aes(fill=lambda_control,color=DSR),size=0.5,stroke=0.1,pch=21),dpi=600) +
    scale_color_manual(values=c('gray','#FC9272')) +
    scale_fill_manual(values=c('gray','#756BB1')) +
    geom_hline(yintercept=c(-1,1),linetype='dashed') +
    scale_x_continuous(trans='log10',limits=c(0.01,15000)) +
    ylim(c(-4,4)) +
    labs(color='lambda control') +
    xlab('mean susceptibility') +
    ylab('log fold-change') +
    guides(color='none',fill='none') +
    theme_bw(base_size=7)
  pv2 <- ggplot(as.data.frame(resASH),aes(x=log2FoldChange,y=minus_log10_pval)) +
    rasterize(geom_point(aes(fill=lambda_control,color=DSR),size=0.5,stroke=0.1,pch=21),dpi=600) +
    scale_color_manual(values=c('gray','#FC9272')) +
    scale_fill_manual(values=c('gray','#756BB1')) +
    xlim(c(-4,4)) +
    ylim(c(0,20)) +
    ylab('-log10(pvalue)') +
    geom_hline(yintercept = -log10(DSR_thresh[1]), linetype='dashed') +
    geom_vline(xintercept = c(-1*DSR_thresh[2],DSR_thresh[2]), linetype='dashed') + 
    theme_bw(base_size=7)
  DSR_plots <- ggarrange(pm1,pv2,ncol=2,widths=c(1,1.5))
  # -------------------------------------------------------------------------- #
  #pdf(file=paste0('figures/DSRs/ma_vol_2_',Contrast[2],ifelse(include_spike_reads,'_spike',''),'.pdf'),width=4,height=1.6)
    annotate_figure(DSR_plots,top=paste(Contrast[2],'vs',Contrast[3]))
  #dev.off()
  # -------------------------------------------------------------------------- #
  
  
  
  # Circos plots with DSR locations
  # Define genome tiles
  dsr <- GRanges(row.names(resASH))
  mcols(dsr) <- DataFrame(resASH)
  sl <- seqlengths(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
  sl <- sl[seqlevels(dsr)]
  sl <- sl[complete.cases(sl)]
  tg <- tileGenome(sl,tilewidth=10e6,cut.last.tile.in.chrom=TRUE)
  #  Count DSRs in genome tiles
  bed <- count_DSRs_in_tiles(dsr=dsr,tiles=tg,padj_thresh=DSR_thresh[1],log2fc_thresh=DSR_thresh[2])
  # Circos plot
  # -------------------------------------------------------------------------- #
  #pdf(file=paste0('figures/DSRs/circos2_',Contrast[2],'.pdf'),width=1.6,height=1.6)
    DSR_plot3(bed=bed,col_range=c(0,30),ylim=c(0,50),outline_points=F,hline_interval=10,cex_modifier=0.5,bg.border=NA)
  #dev.off()
  # -------------------------------------------------------------------------- #
  
  
  
  # Export supplemental files
  #write.csv(resASH,file=paste0(Contrast[2],'-DSR_results.csv'))
  
  
  # Bin Size selection
    library(stringr)
    library(ggplot2)
    se_files <- dir('~/Desktop/trad/figures/DSRs',pattern='se_binSize*',full.names=TRUE)
    pclus_files <- dir('~/Desktop/trad/figures/DSRs',pattern='pclus_*',full.names=TRUE)
    
    se_files <- str_sort(se_files,numeric=TRUE)
    pclus_files <- str_sort(pclus_files,numeric=TRUE)
    
    options(scipen=999)
    rep_cor <- sapply(se_files,function(i) {
      sei <- readRDS(i)
      ci <- assay(sei)
      splt <- strsplit(colnames(ci),'_UV[B|C]_')
      rep <- sapply(splt,'[',2)
      ids <- sapply(splt,'[',1)
      corsi <- sapply(ids[seq(1,length(ids),by=2)],function(j) {
        dati <- ci[,grepl(j,colnames(ci))]
        cor(dati[,1],dati[,2],method='spearman')
      })
      corsi
    },USE.NAMES=TRUE,simplify=FALSE)
    rep_cor <- do.call(rbind,rep_cor)
    
    rep_cor <- as.data.frame.table(rep_cor)
    rep_cor$binSize <- as.numeric(tools::file_path_sans_ext(sapply(strsplit(as.character(rep_cor$Var1),'_binSize_'),'[',2)))
    rep_cor$binSize <- factor(rep_cor$binSize,levels=sort(unique(rep_cor$binSize)))
    
    ggplot(rep_cor,aes(x=binSize,y=Freq)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    scree <- sapply(pclus_files,function(i) {
      pci <- readRDS(i)
      pci$group <- tools::file_path_sans_ext(basename(i))
      pci
    },USE.NAMES=FALSE,simplify=FALSE)
    scree <- as.data.frame(do.call(rbind,scree))
    scree$group <- factor(scree$group,levels=str_sort(unique(scree$group),numeric=TRUE))
    scree$y <- 10^scree$y
    
    ggplot(scree,aes(x=x,y=y,col=group)) +
      geom_line() +
      scale_y_log10() +
      labs(x='Dimensions',y='Percentage of explained variances') +
      scale_color_brewer(palette='Blues') +
      theme_bw()
  