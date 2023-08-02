  library(data.table)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(DESeq2)
  library(pheatmap)
  library(grid)
  library(RColorBrewer)
  library(ggplot2)
  library(tidyr)
    
  
  
  #setwd('/work/data/trad/chip')
  setwd('~/Desktop/tp53/chip/')
  source('~/Documents/GitHub/tradtools/manuscript/old/aggregate_bigwig.R')


  
  # Import top DSRs
  chr <- c(paste0('chr',1:22),'chrX')
  dsr <- read.csv('~/Desktop/trad/DSR_results.csv')
  dsr.gr <- GRanges(dsr[,1])
  mcols(dsr.gr) <- dsr[,-1]
  top_dsr <- dsr.gr[dsr.gr$DSR==TRUE]
  
  
  
  # Define parameters
  dsr_size <- 50e3
  bin_size <- 10e3
  flanking <- 500e3
  bigwig_files <- dir('E017/pval',full.names=TRUE,pattern='.pval.signal.bigwig')
  
  
  
  # Roadmap epigenomics data files
  dat <- sapply(bigwig_files,function(i) {
    aggregate_bigwig(bed_file=top_dsr,bigwig_file=i,bin_size=bin_size,flanking=flanking,ignore.strand=TRUE,verbose=TRUE)
  },simplify=TRUE,USE.NAMES=TRUE)
  colnames(dat) <- sapply(strsplit(colnames(dat),'-|(.pval.signal.bigwig)'),'[',2)
  # saveRDS(dat,'dat_070123.RDS')
  #dat <- readRDS('dat_070123.RDS')
  
  
  
  # Gene density
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  tx <- TxDb.Hsapiens.UCSC.hg19.knownGene
  gene_dens <- coverage(transcripts(tx))
  gene_dens <- gene_dens[names(gene_dens) %in% chr]
  gd <- aggregate_bigwig(top_dsr,bigwig_file=gene_dens,bin_size=bin_size,flanking=flanking)
  # Heatmap
  mat <- t(dat)
  cn <- as.character(1:ncol(mat))
  colnames(mat) <- cn
  COLORS <- rev(RColorBrewer::brewer.pal(9,'RdYlBu'))
  up_num <- round(flanking/bin_size)
  dsr_num <- round(dsr_size/bin_size)
  dn_num <- round(flanking/bin_size)
  col_data <- data.frame(DSR=c(rep('flanking',up_num),c(rep('DSR',dsr_num),rep('flanking',dn_num))),
                         gene_density=gd)
  row.names(col_data) <- cn
  row_data <- data.frame(group=c(rep('A',1),
                                 rep('B',4),
                                 rep('C',1),
                                 rep('D',6),
                                 rep('E',2),
                                 rep('F',1),
                                 rep('G',12),
                                 rep('H',1),
                                 rep('I',1)))
  row.names(row_data) <- c('H3K36me3',
                           'H4K20me1','H3K79me2','H3K79me1','H3K9me1',
                           'DNase',
                           'H3K4me3','H3K9ac','H3K56ac','H2A.Z','H2AK9ac','H2BK5ac',
                           'H3K4me2','H3K18ac',
                           'H3K4me1',
                           'H3K27ac','H4K5ac','H4K8ac','H3K4ac','H3K14ac','H3K23ac','H2AK5ac','H4K91ac','H2BK120ac','H2BK12ac','H2BK15ac','H2BK20ac',
                           'H3K27me3',
                           'H3K9me3')
  
  ph <- pheatmap(mat,cluster_cols = F,scale='row',border_color=NA,breaks=seq(-2,2,length.out=length(COLORS)+1),color=COLORS,fontsize=4,
                 treeheight_row=0,
                 annotation_col=col_data,
                 annotation_row=row_data,
                 annotation_names_row=FALSE,
                 annotation_names_col=FALSE,
                 annotation_colors=list(group=setNames(brewer.pal(9,'Set3'),nm=LETTERS[1:9]),
                                        DSR=setNames(c('gray90','black'),c('flanking','DSR'))),
                 annotation_legend=FALSE,
                 show_rownames=TRUE,
                 show_colnames=FALSE)
  
  pdf('chromatin_heatmap2.pdf',width=4.6,height=1.8)
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
  dev.off()
  
  
  
  # Chromatin states
  # Import data
  cs <- import.bed('~/Desktop/tp53/E017_15_coreMarks_mnemonics.bed.gz')
  seqlevels(cs) <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)
  seqinfo(cs) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
  cs <- keepSeqlevels(cs,value=chr,pruning.mode='coarse')
  cs.names <- unique(cs$name)
  
  dat.cs <- sapply(cs.names,function(i) {
    csi <- cs[cs$name==i]
    csi <- coverage(csi)[chr]
    aggregate_bigwig(bed_file=top_dsr,bigwig_file=csi,bin_size=bin_size,flanking=flanking,ignore.strand=TRUE,verbose=TRUE)
  },simplify=TRUE,USE.NAMES=TRUE)
  
  mat.cs <- t(dat.cs)
  COLORS2 <- rev(RColorBrewer::brewer.pal(9,'PRGn'))
  ph.cs <- pheatmap(mat.cs,cluster_cols = F,scale='row',border_color=NA,breaks=seq(-5,5,length.out=length(COLORS)+1),color=COLORS2,fontsize=4,
                    treeheight_row=0,
                    show_rownames=TRUE,
                    show_colnames=FALSE)
  
  pdf('chromState_heatmap2.pdf',width=4.6,height=0.94)
    grid::grid.newpage()
    grid::grid.draw(ph.cs$gtable)
  dev.off()
  
  
  
  # Dot plot
  dp.dat <- lapply(cs.names,function(i) {
    cs.i <- cs[cs$name==i]
    cs.i <- cs.i[seqnames(cs.i) %in% chr]
    seqlevels(cs.i) <- seqlevelsInUse(cs.i)
    # sanity check
    stopifnot(isDisjoint(cs.i))
    # find fraction of bps annotated to a chromatin state that overlap top DSRs
    total_cs_bases <- sum(width(cs.i))
    cs_bases_overlapping_dsr <- sum(sum(coverage(cs.i)[top_dsr]))
    frac <- cs_bases_overlapping_dsr / total_cs_bases
    # median lesion log2FC in 50kb windows that overlap chromatin state
    ol <- findOverlaps(top_dsr,cs.i,type='any')
    hits <- top_dsr[unique(queryHits(ol))]
    meds <- median(hits$minus_log10_pval)
    list(frac=frac,meds=meds)
  })
  dp.dat2 <- as.data.frame(do.call(rbind,dp.dat))
  dp.dat2 <- as.data.frame(apply(dp.dat2,2,unlist))
  dp.dat2$chromatin_state <- factor(cs.names)
  dp.dat2$chromatin_state <- factor(dp.dat2$chromatin_state,levels=dp.dat2[order(dp.dat2$frac),]$chromatin_state)
  dp.dat2$frac <- as.numeric(dp.dat2$frac)
  dp.dat2$meds <- as.numeric(dp.dat2$meds)
  pdf('dot_plot.pdf',width=3,height=1.4)
  ggplot(dp.dat2,aes(x=chromatin_state,y=frac)) +
    geom_point(aes(color=meds),size=2.5) + 
    scale_color_continuous(type='viridis') + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5))
  dev.off()
  

  
  # Nuclear Lamina
  # Files downloaded from GEO and uniformly re-processed using ENCODE-DCC chip-seq-pipeline2
  # GSE54332: LDF_LMNA:       SRR1141000.srt.nodup_x_SRR1141001.srt.nodup.pval.signal.bigwig
  # GSE81671: HSF_CTL-1_LMNA: SRR3554813.srt.nodup_x_SRR3554814.srt.nodup.pval.signal.bigwig
  # GSE81671: HSF_CTL-2_LMNA: SRR3554823.srt.nodup_x_SRR3554824.srt.nodup.pval.signal.bigwig
  # GSE81671: HSF_CTL-3_LMNA: SRR3554817.srt.nodup_x_SRR3554818.srt.nodup.pval.signal.bigwig
  # GSE49341: LMNB1: SRR944678-SRR944679-SRR944680_x_SRR944690-SRR944691-SRR944692.pval.signal.bigwig
  
  lamin_bigwigs <- c('~/Desktop/lamin/signal/SRR1141000.srt.nodup_x_SRR1141001.srt.nodup.pval.signal.bigwig',
                     '~/Desktop/lamin/signal/SRR3554813.srt.nodup_x_SRR3554814.srt.nodup.pval.signal.bigwig',
                     '~/Desktop/lamin/signal/SRR3554823.srt.nodup_x_SRR3554824.srt.nodup.pval.signal.bigwig',
                     '~/Desktop/lamin/signal/SRR3554817.srt.nodup_x_SRR3554818.srt.nodup.pval.signal.bigwig',
                     '~/Desktop/lamin/signal/SRR944678-SRR944679-SRR944680_x_SRR944690-SRR944691-SRR944692.pval.signal.bigwig')
  lamin_ids <- c('LDF_LMNA',
                 'HSF_CTL-1_LMNA',
                 'HSF_CTL-2_LMNA',
                 'HSF_CTL-3_LMNA',
                 'LMNB1')
    
  lamin <- sapply(lamin_bigwigs,function(i) {
    aggregate_bigwig2(bed_file=top_dsr,bigwig_file=i,bin_size=bin_size,flanking=flanking,ignore.strand=TRUE,verbose=TRUE)
  },simplify=FALSE,USE.NAMES=TRUE)
  names(lamin) <- lamin_ids
  #saveRDS(lamin,'lamin.RDS')
    
    
  # plot_aggregate_bigwig2(lamin$LMNB1,col='purple')
  # COLORS3 <- rev(RColorBrewer::brewer.pal(9,'Spectral'))
  # pheatmap(lamin$LMNB1$matrix[order(rowSums(lamin$LMNB1$matrix),decreasing=TRUE),],
  #          cluster_rows=F,
  #          cluster_cols=F,
  #          scale='none',
  #          legend = F,
  #          color=COLORS3,
  #          breaks=seq(0,50e3,
  #                     length.out=length(COLORS3)+1))
  # save 590 px wide 2950 tall for ~ 300 dpi at 5cm
  
  
  COLORS_DISC <- c('#3182BD','#BAE4B3','#74C476','#31A354','#E6550D')
  pdf(file='LDF_LMNA_profile.pdf',width=2,height=1)
    plot_aggregate_bigwig2(lamin$LDF_LMNA,col=COLORS_DISC[1]) + expand_limits(y=c(4050,4500))
  dev.off()
  plot_aggregate_bigwig2_heatmap(lamin$LDF_LMNA,breaks_type='quantile')
  
  pdf(file='HSF_CTL-1_LMNA_profile.pdf',width=2,height=1)
    plot_aggregate_bigwig2(lamin$`HSF_CTL-1_LMNA`,col=COLORS_DISC[2]) + expand_limits(y=c(3350,3900))
  dev.off()
  plot_aggregate_bigwig2_heatmap(lamin$`HSF_CTL-1_LMNA`,breaks_type='quantile')
  
  pdf(file='HSF_CTL-2_LMNA_profile.pdf',width=2,height=1)
    plot_aggregate_bigwig2(lamin$`HSF_CTL-2_LMNA`,col=COLORS_DISC[3]) + expand_limits(y=c(4100,4500))
  dev.off()
  plot_aggregate_bigwig2_heatmap(lamin$`HSF_CTL-2_LMNA`,breaks_type='quantile')
  
  pdf(file='HSF_CTL-3_LMNA_profile.pdf',width=2,height=1)
    plot_aggregate_bigwig2(lamin$`HSF_CTL-3_LMNA`,col=COLORS_DISC[4]) + expand_limits(y=c(3500,5000))
  dev.off()
  plot_aggregate_bigwig2_heatmap(lamin$`HSF_CTL-3_LMNA`,breaks_type='quantile')
  
  pdf(file='LMNB1_profile.pdf',width=2,height=1)
    plot_aggregate_bigwig2(lamin$LMNB1,col=COLORS_DISC[5]) + expand_limits(y=c(5400,7500))
  dev.off()
  plot_aggregate_bigwig2_heatmap(lamin$LMNB1,breaks_type='sequence',scale_factor=2)
  
  
 
  
  
  

  
  
  
  
  
  # plot_aggregate_bigwig2_heatmap <- function(x) {
  #   CI <- as.data.frame(t(x$CI))
  #   CI$idx <- 1:nrow(CI)
  #   md <- x$metadata
  #   up <- CI$idx[1]
  #   dn <- tail(CI$idx,1)
  #   start <- floor(md$flanking / md$bin_size) + 1
  #   end <- dn - start + 1
  #   mat <- x$matrix
  #   mat <- mat[order(rowSums(mat),decreasing=FALSE),]
  #   df <- as.data.frame.table(mat)
  #   BREAKS <- levels(df$Var2)[c(up,start,end,dn)]
  #   ggplot(df,aes(x=Var2,y=Var1)) +
  #     geom_tile(aes(fill=Freq)) +
  #     scale_x_discrete(breaks=BREAKS,
  #                      labels=c(paste0('-',format_bp(md$flanking)),'s','e',format_bp(md$flanking))) +
  #     scale_fill_gradient2(low='#3288BD',mid='#66C2A5',high='#D53E4F',
  #       limits=c(-1,50e3),oob=scales::squish) +
  #     #scale_fill_manual(breaks=seq(0,50e3,length.out=length(COLORS3)+1)) +
  #     #scale_fill_distiller() +
  #     theme(axis.text.y=element_blank(),
  #           axis.title=element_blank())
  #       
  # 
  #   
  # }
  
    
  
  LDF_LMNA
  
  LDF_LMNA <- import.bed('GSM1313397_HSF_lonzacc2511_LMNA.bed.gz')
  LDF_LMNA <- keepSeqlevels(LDF_LMNA,value=chr,pruning.mode='coarse')
  
  AD04_LMNA.r1 <- import.bed('GSM1313399_HSF_AD04_LMNA_rep1.bed.gz')
  AD04_LMNA.r1 <- keepSeqlevels(AD04_LMNA.r1,value=chr,pruning.mode='coarse')
  
  AD04_LMNA.r2 <- import.bed('GSM1313400_HSF_AD04_LMNA_rep2.bed.gz')
  AD04_LMNA.r2 <- keepSeqlevels(AD04_LMNA.r2,value=chr,pruning.mode='coarse')
  
  

  dat2 <- aggregate_bigwig2(resize(LDF_LMNA,width=1e6,fix='center'), bigwig_file=coverage(dsr.gr,weight='minus_log10_pval'),bin_size=10000,flanking=1e6)
  
  #dat2 <- aggregate_bigwig2(resize(AD04_LMNA.r1,width=1e6,fix='center'), bigwig_file=coverage(dsr.gr,weight='minus_log10_pval'),bin_size=10000,flanking=1e6)
  
  #dat2 <- aggregate_bigwig2(resize(AD04_LMNA.r2,width=1e6,fix='center'), bigwig_file=coverage(dsr.gr,weight='minus_log10_pval'),bin_size=10000,flanking=1e6)
  
  
  #dat2 <- aggregate_bigwig2(top_dsr, bigwig_file=coverage(LDF_LMNA),bin_size=10000,flanking=1e6)
  
  
  COLORS3 <- rev(RColorBrewer::brewer.pal(9,'Spectral'))
  ph.lmna <- pheatmap(dat2[order(rowSums(dat2),decreasing=TRUE),],cluster_rows = F, cluster_cols = F,scale='none',legend = F,color=COLORS3,breaks=seq(0,50e3,length.out=length(COLORS3)+1))
  tiff('LDF_LMNA.tiff',width=0.8,height=2,units='in',res=600)
    grid::grid.newpage()
    grid::grid.draw(ph.lmna$gtable)
  dev.off()
  plot(colSums(dat2),pch=19)
  
  
  pheatmap(dat2,cluster_rows = T, cluster_cols = F,scale='none',treeheight_row = 0)
  
  # composite_plot(LDF_LMNA,var=coverage(dsr.gr,weight='minus_log10_pval'),num_bins=100,upstream=1e6,downstream=1e6)
