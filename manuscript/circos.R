
  library(circlize)
  library(GenomicRanges)
  library(data.table)
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  # Import results from DSR identification
  res <- readRDS('~/Desktop/trad/TP53_resASH.RDS') # see pca.R for details
  dsr <- GRanges(row.names(res))
  mcols(dsr) <- DataFrame(res)
  

  
  # Define genome tiles
  sl <- seqlengths(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
  sl <- sl[seqlevels(dsr)]
  tg <- tileGenome(sl,tilewidth=10e6,cut.last.tile.in.chrom=TRUE)

  
  
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
    circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22),plotType=c('ideogram','labels'))
    # up
    bed_up <- bed[,1:4]
    color_ramp_up <- brewer.pal(9,'Reds')
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
    color_ramp_dn <- brewer.pal(9,'Blues')
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
  
  
  
  #  Count DSRs in genome tiles
  bed <- count_DSRs_in_tiles(dsr=dsr,tiles=tg,padj_thresh=0.05,log2fc_thresh=1)
  
  
  # TP53 col_range=c(0,30); ylim=c(0,50)
  DSR_plot3(bed=bed,col_range=c(0,30),ylim=c(0,50),outline_points=F,hline_interval=5,cex_modifier=0.8,bg.border=NA)
  

  
  
  
  