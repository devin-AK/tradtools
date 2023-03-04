
  library(rtracklayer)
  library(GenomicRanges)
  library(data.table)
  library(pheatmap)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(ggplot2)
  
  
  # Import data
  # chromatin states
  si <- seqinfo(Hsapiens)
  si <- si[seqnames(si)[1:23]]
  cs <- import.bed('~/Desktop/tp53/E017_15_coreMarks_mnemonics.bed.gz')
  cs <- keepSeqlevels(cs,value=seqlevels(si),pruning.mode='coarse')
  seqinfo(cs) <- si
  cs.names <- unique(cs$name)
  #resASH <- read.csv('~/Desktop/trad/DSR_TP53_50kb_resASH.csv')
  resASH <- read.csv('~/Desktop/trad/DSR_100UVC_50kb_resASH.csv')
  resASH <- GRanges(resASH$X,mcols=resASH[,-1])
  resASH <- keepSeqlevels(resASH,value=seqlevels(si),pruning.mode='coarse')
  seqinfo(resASH) <- si
  
  
  # Define top DSRs
  dsr <- resASH[complete.cases(mcols(resASH)),]
  dsr <- dsr[order(dsr$mcols.padj),]
  #top_dsr <- dsr[dsr$padj < 0.0001,]
  dsr.up <- dsr[dsr$mcols.log2FoldChange>0,]
  top_dsr <- dsr[1:10000]
  sl <- seqlevels(si)
  
  
  # Dot Plot
  dp.dat <- lapply(cs.names,function(i) {
    cs.i <- cs[cs$name==i]
    cs.i <- cs.i[seqnames(cs.i) %in% sl]
    seqlevels(cs.i) <- seqlevelsInUse(cs.i)
    seqinfo(cs.i) <- seqinfo(dsr)
    # sanity check
    stopifnot(isDisjoint(cs.i))
    # find fraction of bps annotated to a chromatin state that overlap top DSRs
    total_cs_bases <- sum(width(cs.i))
    cs_bases_overlapping_dsr <- sum(sum(coverage(cs.i)[top_dsr]))
    frac <- cs_bases_overlapping_dsr / total_cs_bases
    # median lesion log2FC in 50kb windows that overlap chromatin state
    ol <- findOverlaps(dsr,cs.i,type='any')
    hits <- dsr[unique(queryHits(ol))]
    meds <- median(hits$mcols.log2FoldChange)
    list(frac=frac,meds=meds)
  })
  dp.dat2 <- as.data.frame(do.call(rbind,dp.dat))
  dp.dat2 <- as.data.frame(apply(dp.dat2,2,unlist))
  dp.dat2$chromatin_state <- factor(cs.names)
  dp.dat2$chromatin_state <- factor(dp.dat2$chromatin_state,levels=dp.dat2[order(dp.dat2$frac),]$chromatin_state)
  dp.dat2$frac <- as.numeric(dp.dat2$frac)
  dp.dat2$meds <- as.numeric(dp.dat2$meds)
  ggplot(dp.dat2,aes(x=chromatin_state,y=frac)) +
    geom_point(aes(color=meds,size=10)) + 
    scale_color_continuous(type='viridis') + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  
  # Heatmap
  dat <- lapply(cs.names,function(i) {
    cs.cov <- coverage(cs[cs$name==i])
    cs.cov <- cs.cov[sl]
    composite_plot(x=top_dsr,var=cs.cov,num_bins=10,upstream=500e3,downstream=500e3,ignore.strand = T)
  })
  names(dat) <- cs.names
  mat <- do.call(rbind,dat)
  pheatmap(mat,cluster_cols=F,scale='row')
  
  
  
  # Profile plot
  composite_plot(x=top_dsr,var=coverage(cs[cs$name=='9_Het'])[sl],num_bins=10,upstream=500e3,downstream=500e3,ignore.strand=T,spar=1)
  
  
  
  

  
  
  
  
  
  
  
  
  