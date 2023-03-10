  
  library(pheatmap)
  library(RColorBrewer)
  bed_file <- '~/Desktop/tp53/TP53_tmp.bed'
  bin_size <- 1000
  flanking <- 50000
  
  # UV damage (TRAD-Seq)
  wt.1 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/trad/results/bw/IMR_200J_UVB_1_possort_markdup_spec1.bw',bin_size=bin_size,flanking=flanking)
  wt.2 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/trad/results/bw/IMR_200J_UVB_2_possort_markdup_spec1.bw',bin_size=bin_size,flanking=flanking)
  tp53.1 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/trad/results/bw/TP53_200J_UVB_1_possort_markdup_spec1.bw',bin_size=bin_size,flanking=flanking)
  tp53.2 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/trad/results/bw/TP53_200J_UVB_2_possort_markdup_spec1.bw',bin_size=bin_size,flanking=flanking)
  kmt5bc.1 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/trad/results/bw/KMT5BC_200J_UVB_1_possort_markdup_spec1.bw',bin_size=bin_size,flanking=flanking)
  kmt5bc.2 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/trad/results/bw/KMT5BC_200J_UVB_2_possort_markdup_spec1.bw',bin_size=bin_size,flanking=flanking)
  wi.1 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/trad/results/bw/WI38_200J_UVB_1_possort_markdup_spec1.bw',bin_size=bin_size,flanking=flanking)
  wi.2 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/trad/results/bw/WI38_200J_UVB_2_possort_markdup_spec1.bw',bin_size=bin_size,flanking=flanking)
  
  
  
  # Chromatin Features (ChIP-Seq)
  wt.laminAC   <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/tp53/chip/Pooled.IMR90LaminACxCtrl.pval.signal.bigwig',bin_size=bin_size,flanking=flanking)
  mut.laminAC  <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/tp53/chip/Pooled.TP53LaminACxCtrl.pval.signal.bigwig',bin_size=bin_size,flanking=flanking)
  wt.laminB    <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/tp53/chip/Pooled.IMRLaminBxCtrl.pvalue.signal.bigwig',bin_size=bin_size,flanking=flanking)
  mut.laminB   <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/tp53/chip/Pooled.TP53LaminBxCtrl.pval.signal.bigwig',bin_size=bin_size,flanking=flanking)
  wt.lbr       <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/tp53/chip/IMR90_LBR_x_ctl.pooled.pval.signal.bigwig',bin_size=bin_size,flanking=flanking)
  mut.lbr      <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/tp53/chip/TP53_LBR_x_ctl.pooled.pval.signal.bigwig',bin_size=bin_size,flanking=flanking)
  wt.h4k20me3  <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/tp53/chip/IMR90H4K20me3_x_ctl.srt.nodup.pval.signal.bigwig',bin_size=bin_size,flanking=flanking)
  mut.h4k20me3 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/tp53/chip/TP53H4K20me3_x_ctl.pooled.pval.signal.bigwig',bin_size=bin_size,flanking=flanking)
  
  save.image(file=paste0('~/Desktop/tp53/kmt_image_binSize',bin_size,'_flaking',as.character(flanking),'.RData'))
  
  
  # Plots
  cpd <- rbind(wt.1,wt.2,tp53.1,tp53.2,kmt5bc.1,kmt5bc.2)
  pheatmap(cpd,cluster_rows=TRUE,cluster_cols=FALSE,scale='row',color=colorRampPalette(rev(brewer.pal(n=7,name='PRGn')))(100),border_color=NA)
  cpd.scale <- t(scale(t(cpd),center=TRUE,scale=TRUE))
  pheatmap(cpd.scale,cluster_rows=TRUE,cluster_cols=FALSE,scale='none',color=colorRampPalette(rev(brewer.pal(n=7,name='PRGn')))(100))
  
  
  tp53diff <- log2((cpd['tp53.1',] + cpd['tp53.2',]) / (cpd['wt.1',] + cpd['wt.2',]))
  kmtdiff  <- log2((cpd['kmt5bc.1',] + cpd['kmt5bc.2',]) / (cpd['wt.1',] + cpd['wt.2',]))
  wi38diff <- log2((kmt5bc.1+kmt5bc.2) / (wi.1+wi.2))
  
  
  chip <- rbind(wt.laminAC,mut.laminAC,wt.laminB,mut.laminB,wt.lbr,mut.lbr,wt.h4k20me3,mut.h4k20me3)
  chip.scale <- t(scale(t(chip),center=TRUE,scale=TRUE))
  pheatmap(chip.scale[1:6,],cluster_rows=FALSE,cluster_cols=FALSE,scale='none')
  pheatmap(chip.scale[7:8,],cluster_rows=FALSE,cluster_cols=FALSE,scale='none')
  
  
  
  aggregate_bigwig <- function(bed_file, bigwig_file, bin_size=10, flanking=0, ignore.strand=TRUE, verbose=TRUE) {
    require(rtracklayer)
    require(GenomicRanges)
    require(data.table)
    if(!ignore.strand) stop('Strand awareness not currently supported')
    stopifnot(file.exists(bed_file))
    stopifnot(file.exists(bigwig_file))
    # Import bed 
    x <- rtracklayer::import.bed(bed_file)
    wx <- width(x)
    mwx <- median(wx)
    if(!all(wx==mwx)) stop('All ranges in bed_file must be same width')
    # Resize bed (flanking sequence)
    if(flanking != 0) {
      nw <- mwx + 2*flanking
      x <- GenomicRanges::resize(x, width=nw, fix='center', use.names=FALSE, ignore.strand=TRUE)
    } else {
      nw <- mwx
    }
    # Process by chromosome
    chr <- seqlevels(x)
    ntile <- round(nw / bin_size)
    dat.chr <- lapply(chr,function(i) {
      if(verbose) message('Processing ',i,'     \r',appendLF=FALSE)
      flush.console()
      xi <- x[seqnames(x)==i]
      vari <- rtracklayer::import(bigwig_file,which=xi,as='RleList')
      vari <- vari[[i]]
      ti <- tile(ranges(xi,use.names=FALSE),n=ntile)
      ti <- unlist(ti)
      v <- Views(vari,ti)
      vs <- viewSums(v,na.rm=TRUE)
      dt <- data.table(matrix(vs,ncol=ntile,byrow=TRUE))
      unname(colSums(dt))
    })
    res <- colSums(matrix(unlist(dat.chr),ncol=ntile,byrow=TRUE))
    if(verbose) message('\nDONE.')
    return(res)
  }
  
  
  
 