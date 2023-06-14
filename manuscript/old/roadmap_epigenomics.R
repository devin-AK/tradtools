  
  
  library(data.table)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(DESeq2)
  
  #setwd('/work/data/trad/chip')
  setwd('~/Desktop/tp53/chip/')
  
  
  format_bp <- function (x, units = 'auto', digits = 1L, ...) {
    known_units <- c('bp','kb','Mb','Gb')
    unit_weights <- c(1,1e3,1e6,1e9)
    units <- match.arg(units, c('auto', known_units))
    if (units == 'auto') {
      test <- log2(x/unit_weights)
      test[test < 0] <- NA_integer_
      idx <- which.min(test)
    } else {
      idx <- match(toupper(units),toupper(known_units))
      if(is.na(idx)) stop('"',units,'" is not a recognized unit. It must be one of ',paste(known_units,collapse=', '))
    }
    unit <- known_units[idx]
    unit_weight <- unit_weights[idx]
    
    paste(format(round(x/unit_weight,digits=digits),scientific=F),unit)
  }
  .composite_data <- function(x, var, ntile, ignore.strand=F, verbose=T) {
    chr <- intersect(names(var),as.character(runValue(seqnames(x))))
    chr.dat <- lapply(chr,function(i) {
      if(verbose) message('Processing ',i,'     \r',appendLF = FALSE)
      flush.console()
      vari <- var[[i]]
      xi <- x[seqnames(x)==i]
      strnd <- strand(xi)
      tiles <- tile(ranges(xi),n=ntile)
      v <- Views(vari,unlist(unname(tiles)))
      vs <- viewSums(v,na.rm=T)
      dt <- data.table(matrix(vs,ncol=ntile,byrow=TRUE))
      if(!ignore.strand) {
        idx <- as.logical(strnd == '-')
        dt[idx,] <- dt[idx,rev(.SD)]
      }
      #unname(colSums(dt))
      unname(colMeans(dt))
    })
    #res <- colSums(matrix(unlist(chr.dat),ncol=ntile,byrow=TRUE))
    res <- colMeans(matrix(unlist(chr.dat),ncol=ntile,byrow=TRUE))
    if(verbose) message('\nDONE.')
    return(res)
  }
  
  composite_plot <- function(x, var, num_bins=NULL, bin_size=NULL, upstream=0, downstream=0, ignore.strand=F, spar=NULL, verbose=TRUE, xlab='', ylab='', ...) {
    ### checks
    if(!inherits(x,'GRanges')) stop('"x" must be a GRanges object')
    x <- granges(x,use.names=F)
    if(!inherits(var,'RleList')) stop('"var" must be an RleList object, i.e. made using coverage()')
    stopifnot(isTRUEorFALSE(ignore.strand))
    stopifnot(isTRUEorFALSE(verbose))
    if(!all(names(var) %in% seqlevels(x))) stop('All of the names in the coverage vector "var" must be in seqlevels(x)')
    if(!is.null(num_bins) & !is.null(bin_size)) stop('Only one of "num_bins" and "bin_size" can be specified (not NULL)')
    if(all(is.null(num_bins),is.null(bin_size))) {
      if(verbose) message('Using "bin_size" of 100.')
      bin_size <- 100
    }
    if(!is.null(bin_size)) {
      wx1 <- width(x[1])
      if(!all(width(x)==wx1)) stop('If specifying "bin_size", all elements in x must have the same width.')
      ntile <- floor(wx1 / bin_size)
    } else {
      ntile <- num_bins
    }
    if(min(width(x)) < ntile) stop('One or more ranges in "x" has width < "bin_size"')
    ### get the regions flanking x
    up <- trim(GenomicRanges::flank(x,width=upstream,start=T,both=F,ignore.strand=ignore.strand))
    if(ignore.strand) strand(up) <- '+'
    ntile.up <- floor(upstream/(mean(width(x))/ntile))
    dn <- trim(GenomicRanges::flank(x,width=downstream,start=F,both=F,ignore.strand=ignore.strand))
    if(ignore.strand) strand(dn) <- '+'
    ntile.dn <- floor(downstream/(mean(width(x))/ntile))
    idx <- width(x) > ntile & width(up) > ntile.up & width(dn) > ntile.dn
    if(verbose) message('Removed ',sum(!idx),' ranges')
    ### Process
    upstr <-     .composite_data(x=up[idx],var=var,ntile=ntile.up,ignore.strand=ignore.strand,verbose=verbose)
    center <- .composite_data(x=x[idx], var=var, ntile=ntile, ignore.strand=ignore.strand, verbose=verbose)
    dnstr <-     .composite_data(x=dn[idx], var=var,ntile=ntile.dn,ignore.strand=ignore.strand,verbose=verbose)
    ### Plot
    res <- c(upstr,center,dnstr)
    at <- c(1,length(upstr),length(upstr)+length(center),length(res))
    lab <- c(paste0('-',format_bp(upstream)),'start','end',paste0('+',format_bp(downstream)))
    #plot(res,type='l',xaxt='n',ylab='',xlab='',col=rgb(47,79,79,alpha = 125,maxColorValue = 255),...)
    plot(res,type='l',xaxt='n',ylab=ylab,xlab=xlab,col='skyblue1',lwd=0.25,...)
    if(!is.null(spar)) lines(smooth.spline(res,spar=spar),col='tan1',lwd=2)
    axis(1,at=at,labels=lab,las=2)
    abline(v=c(length(upstr),length(upstr)+length(center)),lty=2)
    
    return(res)
    
  }
  
  #dsr <- readRDS('/work/data/trad/TP53_resASH_DSR_50kb.RDS')
  dsr <- readRDS('~/Desktop/trad/TP53_resASH_DSR_50kb.RDS')
  top_dsr <- dsr[!is.na(dsr$padj) & dsr$padj < 0.00001,]
  top_dsr <- as(row.names(top_dsr),'GRanges')
  chr <- c(paste0('chr',1:22),'chrX')
  

  sample <- 'TP53'
  details <- 'UP_DSR_50kb'
  fls <- dir('E017/pval',full.names=T,pattern='.pval.signal.bigwig')
  #fls <- dir(getwd(),full.names=T,pattern='.signal.bigwig')
  fls.bn <- paste(sample,details,sapply(strsplit(basename(fls),'.pval.signal.b'),'[',1),sep='_')
  pdf(file=paste0('plots/',sample,'_',details,'.pdf'))
  dat2 <- lapply(fls,function(i) {
    bn <- sapply(strsplit(basename(i),'.pval.signal.b'),'[',1)
    vari <- rtracklayer::import.bw(i,as='RleList')
    vari <- vari[names(vari) %in% chr]
    final_nm <- paste0(paste(sample,details,bn,length(top_dsr),sep='_'),'.pdf')
    composite_plot(x=top_dsr, var=vari, num_bins=10, upstream=500e3, downstream=500e3, ignore.strand=F, spar=0.95, main=final_nm)
  })
  dev.off()
  dsr_char2 <- do.call('rbind',dat2)
  row.names(dsr_char2) <- fls.bn
  
  
  
  
  
  
  
  
  hs <- BSgenome.Hsapiens.UCSC.hg19
  si <- seqinfo(hs)
  
  WHICH <- as(si,'GRanges')[1:5]
  sample <- 'TP53'
  details <- 'UP_DSR_50kb_chr1-5'
  fls <- dir('E017/pval',full.names=T,pattern='.pval.signal.bigwig')
  fls <- c(fls,dir(getwd(),full.names=T,pattern='.signal.bigwig'))
  chr <- c(paste0('chr',1:22),'chrX')
  fls.bn <- paste(sample,details,sapply(strsplit(basename(fls),'signal.b'),'[',1),sep='_')
  pdf(file=paste0('plots/',sample,'_',details,'.pdf'))
  dat2 <- lapply(fls,function(i) {
    bn <- sapply(strsplit(basename(i),'.pval.signal.b'),'[',1)
    vari <- rtracklayer::import.bw(i,which=WHICH,as='RleList')
    vari <- vari[names(vari) %in% chr]
    final_nm <- paste0(paste(sample,details,bn,length(top_dsr),sep='_'),'.pdf')
    composite_plot(x=top_dsr, var=vari, num_bins=10, upstream=500e3, downstream=500e3, ignore.strand=F, spar=0.95, main=final_nm)
  })
  dev.off()
  dsr_char2 <- do.call('rbind',dat2)
  row.names(dsr_char2) <- fls.bn
  
  
  pheatmap(dsr_char2,cluster_rows = T,cluster_cols = F,scale='row')
  
  tmp <- log2(dsr_char2["TP53_UP_DSR_50kb_chr1-5_TP53H4K20me3_x_ctl.pooled.pval.",] / dsr_char2["TP53_UP_DSR_50kb_chr1-5_IMR90H4K20me3_x_ctl.srt.nodup.pval.",])
  
  plot(tmp)
  
  
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  tx <- TxDb.Hsapiens.UCSC.hg19.knownGene
  gene_dens <- coverage(transcripts(tx))
  gene_dens <- gene_dens[names(gene_dens) %in% chr]
  
  composite_plot(x=top_dsr,var=gene_dens,num_bins=10,upstream=500e3,downstream=500e3,ignore.strand=T,spar=1)
  
  
  kmt1 <- import.bw('~/Desktop/trad/results/bw/KMT5BC_200J_UVB_1_possort_markdup_spec1.bw',as='RleList')
  kmt1 <- kmt1[names(kmt1) %in% chr]
  kmt1dat <- composite_plot(x=top_dsr,var=kmt1,num_bins=10,upstream=500e3,downstream=500e3,ignore.strand=T,spar=1)
  
  kmt2 <- import.bw('~/Desktop/trad/results/bw/KMT5BC_200J_UVB_2_possort_markdup_spec1.bw',as='RleList')
  kmt2 <- kmt2[names(kmt2) %in% chr]
  kmt2dat <- composite_plot(x=top_dsr,var=kmt2,num_bins=10,upstream=500e3,downstream=500e3,ignore.strand=T,spar=1)
  
  wt1 <- import.bw('~/Desktop/trad/results/bw/IMR_200J_UVB_1_possort_markdup_spec1.bw',as='RleList')
  wt1 <- wt1[names(wt1) %in% chr]
  wt1dat <- composite_plot(x=top_dsr,var=wt1,num_bins=10,upstream=500e3,downstream=500e3,ignore.strand=T,spar=1)
  
  wt2 <- import.bw('~/Desktop/trad/results/bw/IMR_200J_UVB_2_possort_markdup_spec1.bw',as='RleList')
  wt2 <- wt2[names(wt2) %in% chr]
  wt2dat <- composite_plot(x=top_dsr,var=wt2,num_bins=10,upstream=500e3,downstream=500e3,ignore.strand=T,spar=1)
  
  wi1 <- import.bw('~/Desktop/trad/results/bw/WI38_200J_UVB_1_possort_markdup_spec1.bw',as='RleList')
  wi1 <- wi1[names(wi1) %in% chr]
  wi1dat <- composite_plot(x=top_dsr,var=wi1,num_bins=10,upstream=500e3,downstream=500e3,ignore.strand=T,spar=1)
  
  wi2 <- import.bw('~/Desktop/trad/results/bw/WI38_200J_UVB_2_possort_markdup_spec1.bw',as='RleList')
  wi2 <- wi2[names(wi2) %in% chr]
  wi2dat <- composite_plot(x=top_dsr,var=wi2,num_bins=10,upstream=500e3,downstream=500e3,ignore.strand=T,spar=1)
  
  
  alldat <- cbind(wt1dat,wt2dat,kmt1dat,kmt2dat)
  pheatmap(cor(alldat))
  
  pheatmap(t(alldat),cluster_rows = T,cluster_cols = F,scale='row')
  
  
  
  
  wi1 <- import.bw('~/Desktop/trad/results/bw/IMR_200J_UVB_2_possort_markdup_spec1.bw',as='RleList')
  wi1 <- wi1[names(wi1) %in% chr]
  wi1dat <- composite_plot(x=top_dsr,var=wi1,num_bins=100,upstream=5000e3,downstream=5000e3,ignore.strand=T,spar=NULL)
  
  
  
  wi1 <- import.bw('~/Desktop/trad/results/bw/KMT5BC_200J_UVB_1_possort_markdup_spec1.bw',as='RleList')
  wi1 <- wi1[names(wi1) %in% chr]
  wi1dat <- composite_plot(x=top_dsr,var=wi1,num_bins=50,upstream=100e3,downstream=100e3,ignore.strand=T,spar=NULL)
  
  
  
  WHICH <- resize(top_dsr,width=200e3,fix='center')
  WHICH <- reduce(WHICH)
  
  wi1 <- import.bw('~/Desktop/tp53/chip/IMR90H4K20me3_x_ctl.srt.nodup.pval.signal.bigwig',which=WHICH,as='RleList')
  wi1 <- wi1[names(wi1) %in% chr]
  wi1dat <- composite_plot(x=top_dsr,var=wi1,num_bins=50,upstream=50e3,downstream=50e3,ignore.strand=T,spar=NULL)
  
  
  wi1dat <- composite_plot(x=top_dsr,var=wi1,num_bins=100,upstream=50e3,downstream=50e3,ignore.strand=T,spar=NULL)
  
  
  
  
  