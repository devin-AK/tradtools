  
  # library(pheatmap)
  # library(RColorBrewer)
  # bed_file <- '~/Desktop/tp53/TP53_tmp.bed'
  # bin_size <- 1000
  # flanking <- 50000
  # 
  # # UV damage (TRAD-Seq)
  # wt.1 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/trad/results/bw/IMR_200J_UVB_1_possort_markdup_spec1.bw',bin_size=bin_size,flanking=flanking)
  # wt.2 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/trad/results/bw/IMR_200J_UVB_2_possort_markdup_spec1.bw',bin_size=bin_size,flanking=flanking)
  # tp53.1 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/trad/results/bw/TP53_200J_UVB_1_possort_markdup_spec1.bw',bin_size=bin_size,flanking=flanking)
  # tp53.2 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/trad/results/bw/TP53_200J_UVB_2_possort_markdup_spec1.bw',bin_size=bin_size,flanking=flanking)
  # kmt5bc.1 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/trad/results/bw/KMT5BC_200J_UVB_1_possort_markdup_spec1.bw',bin_size=bin_size,flanking=flanking)
  # kmt5bc.2 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/trad/results/bw/KMT5BC_200J_UVB_2_possort_markdup_spec1.bw',bin_size=bin_size,flanking=flanking)
  # wi.1 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/trad/results/bw/WI38_200J_UVB_1_possort_markdup_spec1.bw',bin_size=bin_size,flanking=flanking)
  # wi.2 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/trad/results/bw/WI38_200J_UVB_2_possort_markdup_spec1.bw',bin_size=bin_size,flanking=flanking)
  # 
  # 
  # 
  # # Chromatin Features (ChIP-Seq)
  # wt.laminAC   <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/tp53/chip/Pooled.IMR90LaminACxCtrl.pval.signal.bigwig',bin_size=bin_size,flanking=flanking)
  # mut.laminAC  <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/tp53/chip/Pooled.TP53LaminACxCtrl.pval.signal.bigwig',bin_size=bin_size,flanking=flanking)
  # wt.laminB    <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/tp53/chip/Pooled.IMRLaminBxCtrl.pvalue.signal.bigwig',bin_size=bin_size,flanking=flanking)
  # mut.laminB   <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/tp53/chip/Pooled.TP53LaminBxCtrl.pval.signal.bigwig',bin_size=bin_size,flanking=flanking)
  # wt.lbr       <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/tp53/chip/IMR90_LBR_x_ctl.pooled.pval.signal.bigwig',bin_size=bin_size,flanking=flanking)
  # mut.lbr      <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/tp53/chip/TP53_LBR_x_ctl.pooled.pval.signal.bigwig',bin_size=bin_size,flanking=flanking)
  # wt.h4k20me3  <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/tp53/chip/IMR90H4K20me3_x_ctl.srt.nodup.pval.signal.bigwig',bin_size=bin_size,flanking=flanking)
  # mut.h4k20me3 <- aggregate_bigwig(bed_file=bed_file,bigwig_file='~/Desktop/tp53/chip/TP53H4K20me3_x_ctl.pooled.pval.signal.bigwig',bin_size=bin_size,flanking=flanking)
  # 
  # save.image(file=paste0('~/Desktop/tp53/kmt_image_binSize',bin_size,'_flaking',as.character(flanking),'.RData'))
  # 
  # 
  # # Plots
  # cpd <- rbind(wt.1,wt.2,tp53.1,tp53.2,kmt5bc.1,kmt5bc.2)
  # pheatmap(cpd,cluster_rows=TRUE,cluster_cols=FALSE,scale='row',color=colorRampPalette(rev(brewer.pal(n=7,name='PRGn')))(100),border_color=NA)
  # cpd.scale <- t(scale(t(cpd),center=TRUE,scale=TRUE))
  # pheatmap(cpd.scale,cluster_rows=TRUE,cluster_cols=FALSE,scale='none',color=colorRampPalette(rev(brewer.pal(n=7,name='PRGn')))(100))
  # 
  # 
  # tp53diff <- log2((cpd['tp53.1',] + cpd['tp53.2',]) / (cpd['wt.1',] + cpd['wt.2',]))
  # kmtdiff  <- log2((cpd['kmt5bc.1',] + cpd['kmt5bc.2',]) / (cpd['wt.1',] + cpd['wt.2',]))
  # wi38diff <- log2((kmt5bc.1+kmt5bc.2) / (wi.1+wi.2))
  # 
  # 
  # chip <- rbind(wt.laminAC,mut.laminAC,wt.laminB,mut.laminB,wt.lbr,mut.lbr,wt.h4k20me3,mut.h4k20me3)
  # chip.scale <- t(scale(t(chip),center=TRUE,scale=TRUE))
  # pheatmap(chip.scale[1:6,],cluster_rows=FALSE,cluster_cols=FALSE,scale='none')
  # pheatmap(chip.scale[7:8,],cluster_rows=FALSE,cluster_cols=FALSE,scale='none')
  # 
  # 
  # 
  # aggregate_bigwig1 <- function(bed_file, bigwig_file, bin_size=10, flanking=0, ignore.strand=TRUE, verbose=TRUE) {
  #   require(rtracklayer)
  #   require(GenomicRanges)
  #   require(data.table)
  #   if(!ignore.strand) stop('Strand awareness not currently supported')
  #   stopifnot(file.exists(bed_file))
  #   stopifnot(file.exists(bigwig_file))
  #   # Import bed 
  #   x <- rtracklayer::import.bed(bed_file)
  #   wx <- width(x)
  #   mwx <- median(wx)
  #   if(!all(wx==mwx)) stop('All ranges in bed_file must be same width')
  #   # Resize bed (flanking sequence)
  #   if(flanking != 0) {
  #     nw <- mwx + 2*flanking
  #     x <- GenomicRanges::resize(x, width=nw, fix='center', use.names=FALSE, ignore.strand=TRUE)
  #   } else {
  #     nw <- mwx
  #   }
  #   # Process by chromosome
  #   chr <- seqlevels(x)
  #   ntile <- round(nw / bin_size)
  #   dat.chr <- lapply(chr,function(i) {
  #     if(verbose) message('Processing ',i,'     \r',appendLF=FALSE)
  #     flush.console()
  #     xi <- x[seqnames(x)==i]
  #     vari <- rtracklayer::import(bigwig_file,which=xi,as='RleList')
  #     vari <- vari[[i]]
  #     ti <- tile(ranges(xi,use.names=FALSE),n=ntile)
  #     ti <- unlist(ti)
  #     v <- Views(vari,ti)
  #     vs <- viewSums(v,na.rm=TRUE)
  #     dt <- data.table(matrix(vs,ncol=ntile,byrow=TRUE))
  #     unname(colSums(dt))
  #   })
  #   res <- colSums(matrix(unlist(dat.chr),ncol=ntile,byrow=TRUE))
  #   if(verbose) message('\nDONE.')
  #   return(res)
  # }
  
  
  
  aggregate_bigwig <- function(bed_file, bigwig_file, bin_size=10, flanking=0, ignore.strand=TRUE, verbose=TRUE) {
    require(rtracklayer)
    require(GenomicRanges)
    require(data.table)
    # check inputs
    if(!ignore.strand) stop('Strand awareness not currently supported')
    stopifnot('bin_size must be numeric'=is.numeric(bin_size))
    stopifnot('flanking must be numeric'=is.numeric(flanking))
    if(is.character(bigwig_file)) {
      stopifnot('bigwig_file does not exist. Double check PATH and make sure file is present at that location'=file.exists(bigwig_file))
      stopifnot('bigwig_file must be a BigWig file (file extension .bw or .bigwig)'=tolower(tools::file_ext(bigwig_file)) %in% c('bw','bigwig'))
      bigwig_preloaded <- FALSE
    } else {
      stopifnot('If bigwig_file is an R object, it must be RleList object (e.g. made with coverage()'=inherits(bigwig_file,'RleList'))
      bigwig_preloaded <- TRUE
    }
    if(is.character(bed_file)) {
      # Import bed 
      if(verbose) message('Importing BED file ',basename(bed_file))
      x <- rtracklayer::import.bed(bed_file)
    } else {
      x <- bed_file
    }
    stopifnot('supplied bed_file must be GRanges object or a PATH to a valid BED file that can be imported'=inherits(x,'GRanges'))
    # Process
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
      if(bigwig_preloaded) {
        vari <- bigwig_file[[i]]
      } else {
        vari <- rtracklayer::import(bigwig_file,which=xi,as='RleList')
        vari <- vari[[i]]
      }
      ti <- tile(ranges(xi,use.names=FALSE),n=ntile)
      ti <- unlist(ti)
      v <- Views(vari,ti)
      vs <- viewSums(v,na.rm=TRUE)
      dt <- data.table(matrix(vs,ncol=ntile,byrow=TRUE))
      unname(colSums(dt))
    })
    sums <- colSums(matrix(unlist(dat.chr),ncol=ntile,byrow=TRUE))
    avg_per_bin <- sums / length(x)
    # calculate average per bp
    avg_per_bp <- avg_per_bin / bin_size
    if(verbose) message('\nDONE.')
    return(avg_per_bp)
  }
  
  
  
  aggregate_bigwig2 <- function(bed_file, bigwig_file, bin_size=10, flanking=0, ignore.strand=TRUE, verbose=TRUE) {
    require(rtracklayer)
    require(GenomicRanges)
    require(data.table)
    # check inputs
    if(!ignore.strand) stop('Strand awareness not currently supported')
    stopifnot('bin_size must be numeric'=is.numeric(bin_size))
    stopifnot('flanking must be numeric'=is.numeric(flanking))
    if(is.character(bigwig_file)) {
      stopifnot('bigwig_file does not exist. Double check PATH and make sure file is present at that location'=file.exists(bigwig_file))
      stopifnot('bigwig_file must be a BigWig file (file extension .bw or .bigwig)'=tolower(tools::file_ext(bigwig_file)) %in% c('bw','bigwig'))
      bigwig_preloaded <- FALSE
    } else {
      stopifnot('If bigwig_file is an R object, it must be RleList object (e.g. made with coverage()'=inherits(bigwig_file,'RleList'))
      bigwig_preloaded <- TRUE
    }
    if(is.character(bed_file)) {
      # Import bed 
      if(verbose) message('Importing BED file ',basename(bed_file))
      x <- rtracklayer::import.bed(bed_file)
    } else {
      x <- bed_file
    }
    stopifnot('supplied bed_file must be GRanges object or a PATH to a valid BED file that can be imported'=inherits(x,'GRanges'))
    # Process
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
    dat.chr <- sapply(chr,function(i) {
      if(verbose) message('Processing ',i,'     \r',appendLF=FALSE)
      flush.console()
      xi <- x[seqnames(x)==i]
      if(bigwig_preloaded) {
        vari <- bigwig_file[[i]]
      } else {
        vari <- rtracklayer::import(bigwig_file,which=xi,as='RleList')
        vari <- vari[[i]]
      }
      ti <- tile(ranges(xi,use.names=FALSE),n=ntile)
      ti <- unlist(ti)
      v <- Views(vari,ti)
      vs <- viewSums(v,na.rm=TRUE)
      matrix(vs,ncol=ntile,byrow=TRUE)
    })
    mat <- do.call(rbind,dat.chr)
    # Calculate 95% CI
    if(verbose) message('\nCalculating 95% confidence interval')
    .calculate_matrix_columns_CI <- function(x) {
      # 95% CI
      stopifnot(inherits(x,'matrix'))
      apply(x,2,function(x) {
        mu <- mean(x)
        ci <- mu+c(-1.96,1.96)*sd(x)/sqrt(length(x))
        c('mean'=mu,'lower_CI'=ci[1],'upper_CI'=ci[2])
      },simplify=TRUE)
    }
    CI <- .calculate_matrix_columns_CI(mat)
    if(verbose) message('DONE.')
    return(list('matrix'=mat,'CI'=CI,'metadata'=data.frame('bin_size'=bin_size,'flanking'=flanking)))
  }
  
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
  
  plot_aggregate_bigwig2 <- function(x,col='darkred',CI_color='gray80') {
    CI <- as.data.frame(t(x$CI))
    CI$idx <- 1:nrow(CI)
    md <- x$metadata
    up <- CI$idx[1]
    dn <- tail(CI$idx,1)
    start <- floor(md$flanking / md$bin_size) + 1
    end <- dn - start + 1
    ggplot(CI) +
      geom_ribbon(aes(x=idx,ymin=lower_CI,ymax=upper_CI),fill=CI_color) +
      geom_line(aes(x=idx,y=mean),color=col) +
      ylab('Mean') +
      theme_bw() +
      geom_vline(xintercept=c(start,end),linetype='dashed') +
      scale_x_continuous(breaks=c(up,start,end,dn),
                         labels=c(paste0('-',format_bp(md$flanking)),'s','e',format_bp(md$flanking))) +
      theme(axis.title.x=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank())
  }
  
  plot_aggregate_bigwig2_heatmap <- function(x,palette='Spectral',scale=TRUE,breaks_type='sequence',scale_factor=1.2,rev=TRUE) {
    mat <- x$matrix
    if(scale) mat <- scale(mat,center=FALSE)
    .quantile_breaks <- function(xs, n = 10) {
      breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
      breaks[!duplicated(breaks)]
    }
    N=8
    if(breaks_type=='quantile') mat_breaks <- .quantile_breaks(mat, n=N)*scale_factor
    if(breaks_type=='sequence') mat_breaks <- seq(0,quantile(mat,0.97)*scale_factor,length.out=N)
    COLORS <- brewer.pal(N-1,palette)
    if(rev) COLORS <- rev(COLORS)
    pheatmap(mat[order(rowSums(mat),decreasing=TRUE),],
             cluster_rows=F,
             cluster_cols=F,
             scale='none',
             legend = F,
             color=COLORS,
             breaks=mat_breaks)
  }
  
  
  #bed_file <- '~/Desktop/tp53/TP53_tmp.bed'
  #dat <- aggregate_bigwig2(bed_file=bed_file,bigwig_file='~/Desktop/trad/results/bw/IMR_200J_UVB_1_possort_markdup_spec1.bw',bin_size=1000,flanking=50000)
  # 
  
  
  

    
    
  
  