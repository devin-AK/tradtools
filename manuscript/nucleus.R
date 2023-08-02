
  setwd('~/Desktop/tp53/chrom3d/')
  
  
  # Load libraries
  library(GenomicRanges)
  library(rgl)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(viridis)
  library(scales)
  library(RColorBrewer)
  library(ggplot2)
  
  
  
  # Functions
  force_seqinfo <- function(x, reference=NULL, chromosomes_to_keep=NULL) {
    require('GenomicRanges')
    require('GenomeInfoDb')
    stopifnot(is(x,'GRanges'))
    if(is.null(reference)) {
      si <- seqinfo(x)
      ss <- seqlevelsStyle(x)
      sl <- seqlevels(x)
    } else {
      si <- seqinfo(reference)
      ss <- seqlevelsStyle(si)
      seqlevelsStyle(x) <- ss
      sl <- seqlevels(si)
      sl <- intersect(sl,seqlevels(x))
    }
    if(!is.null(chromosomes_to_keep)) {
      sl <- intersect(sl,chromosomes_to_keep)
    }
    seqlevels(si) <- sl
    seqlevels(x,pruning.mode='coarse') <- sl
    seqinfo(x) <- si
    return(x)
  }
  
  import_DSR_stats <- function(x,reference=NULL,chromosomes_to_keep=NULL,na.rm=FALSE) {
    require('GenomicRanges')
    dsr <- read.csv(x,stringsAsFactors=FALSE)
    MCOLS <- dsr[,-1]
    dsr <- GenomicRanges::GRanges(dsr[,1])
    mcols(dsr) <- MCOLS
    dsr <- force_seqinfo(dsr,reference=reference,chromosomes_to_keep=chromosomes_to_keep)
    if(na.rm) {
      dsr <- dsr[complete.cases(mcols(dsr))]
    }
    dsr
  }
  
  process_cmm_file <- function(cmm_file, numeric_cols=c('x','y','z','radius','r','g','b'), genome) {
    require(XML)
    require(rgl)
    require(GenomicRanges)
    xml <- xmlInternalTreeParse(cmm_file)
    dat <- xmlToList(xml)
    nm <- names(table(names(dat))) # children names (e.g. link and marker in cmm file)
    # format data
    dfl <- lapply(nm, function(i) {
      df <- dat[names(dat) == i]
      df <- as.data.frame(do.call('rbind',df),stringsAsFactors = FALSE)
      row.names(df) <- NULL
      df_numeric <- colnames(df) %in% numeric_cols
      df[, df_numeric] <- sapply(df[, df_numeric],as.numeric)
      df
    })
    names(dfl) <- nm
    # prepare colors
    marker_col <- rgb(red=dfl$marker$r, green=dfl$marker$g, blue=dfl$marker$b, maxColorValue=1)
    link_col   <- rgb(red=dfl$link$r, green=dfl$link$g, blue=dfl$link$b, maxColorValue=1)
    # prepare granges
    bead_id <- dfl$marker$beadID
    bead_id2 <- gsub('_[A|B]','',bead_id)
    gr <- GRanges(bead_id2)
    start(gr) <- start(gr)+1
    gr <- sortSeqlevels(gr)
    sl <- seqlevelsInUse(gr)
    seqinfo(gr) <- seqinfo(genome)
    gr <- keepSeqlevels(gr,sl,'coarse')
    res <- list(marker=dfl$marker, link=dfl$link, marker_col=marker_col, link_col=link_col, ranges=gr)
    return(res)
  }
  
  sphere1.f <- function(x0 = 0, y0 = 0, z0 = 0, r = 1, n = 101, ...){
    require('rgl')
    f <- function(s,t){ 
      cbind(   r * cos(t)*cos(s) + x0,
               r *        sin(s) + y0,
               r * sin(t)*cos(s) + z0)
    }
    persp3d(f, slim = c(-pi/2,pi/2), tlim = c(0, 2*pi), n = n, add = T, ...)
  }
  
  
  
  # Define necessary files and parameters
  trad_signal_wt.1 <- '~/Desktop/trad/results/bw/IMR_200J_UVB_1_possort_markdup_spec1.bw' # Generated in this study
  trad_signal_wt.2 <- '~/Desktop/trad/results/bw/IMR_200J_UVB_2_possort_markdup_spec1.bw' # Generated in this study
  trad_signal_mut.1 <- '~/Desktop/trad/results/bw/TP53_200J_UVB_1_possort_markdup_spec1.bw' # Generated in this study
  trad_signal_mut.2 <- '~/Desktop/trad/results/bw/TP53_200J_UVB_2_possort_markdup_spec1.bw' # Generated in this study
  cmm_file <- '~/Desktop/tp53/chrom3d/IMR90_inter_intra_chr_w_LADs.diploid.cmm'
  REF <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
 
  
  
  # Process
  get_sum_in_ranges <- function(gr,bigwig_file,reference,chromosomes_to_keep) {
    bw <- import.bw(bigwig_file)
    bw <- force_seqinfo(bw,reference=reference,chromosomes_to_keep=chromosomes_to_keep)
    bw <- coverage(bw,weight='score')
    width(gr) * binnedAverage(bins=gr,numvar=bw,varname='x')$x
  }
  cmm <- process_cmm_file(cmm_file,genome=REF)
  chr_to_keep <- seqlevels(cmm$ranges)
  d1 <- get_sum_in_ranges(cmm$ranges,bigwig_file=trad_signal_wt.1,reference=REF,chromosomes_to_keep=chr_to_keep)
  d2 <- get_sum_in_ranges(cmm$ranges,bigwig_file=trad_signal_wt.2,reference=REF,chromosomes_to_keep=chr_to_keep)
  d3 <- get_sum_in_ranges(cmm$ranges,bigwig_file=trad_signal_mut.1,reference=REF,chromosomes_to_keep=chr_to_keep)
  d4 <- get_sum_in_ranges(cmm$ranges,bigwig_file=trad_signal_mut.2,reference=REF,chromosomes_to_keep=chr_to_keep)
  x <- log2((d3+d4)/(d1+d2))
  ncolors <- 5
  pal <- viridis(ncolors,option='D')
  #pal <- rev(brewer.pal(ncolors,'RdBu'))
  Breaks <- c(-Inf,seq(from=quantile(x,na.rm=T,0.25),to=quantile(x,na.rm=T,0.75),length.out=ncolors-1),Inf)
  cuts <- cut(x,breaks=Breaks)
  Cols <- pal[cuts]
  dat <- cbind(cmm$marker[,c('x','y','z','radius','beadID')],colors=Cols,lfc=x)
  dat$d <- sqrt((dat$x)^2 + (dat$y)^2 + (dat$z)^2)
  dat <- dat[complete.cases(dat),]
  dat <- dat[is.finite(dat$lfc),]
  
  
  
  ### Plot
  # association
  pdf('~/Desktop/radial.pdf',width=1.5,height=1.5)
  ggplot(dat,aes(x=d,y=lfc)) +
    #geom_point() +
    geom_smooth(method='loess') +
    xlab(expression(paste('Radial Position (',mu,'m)'))) +
    ylab('Differential Susc. (log2-FC)') +
    theme_bw()
  dev.off()
  
  # legend
  show_col(pal,ncol=ncolors,label=FALSE)
  
  # 3d model
  open3d()
  chromatin <- plot3d(x=dat$x,y=dat$y,z=dat$z,type='s',radius=dat$radius,col=dat$colors,
                     axes=FALSE,
                     box=FALSE,
                     xlab='',
                     ylab='',
                     zlab='',
                     xlim=c(-6.5,6.5),ylim=c(-6.5,6.5),zlim=c(-6.5,6.5))
  cross_section <- planes3d(a=0,b=1,c=0,d=0,col='black',alpha=0.8)
  nuclear_envelope <- sphere1.f(r=6.2,col='tan')
  root <- currentSubscene3d()
  newSubscene3d('inherit','inherit','inherit',copyShapes=TRUE,parent=root)
  clipplanes3d(0,1,0,0)
  useSubscene3d(root)
  delFromSubscene3d(nuclear_envelope)
  um <- matrix(data=c(-0.8,-0.5,-0.3,0,0.6,-0.7,-0.2,0,-0.1,-0.3,1,0,0,0,0,1),byrow=TRUE,ncol=4)
  view3d(userMatrix=um,zoom=0.5)
  snapshot3d(filename='chrom3d_model.png',fmt='png',width=500,height=500,webshot=FALSE)
  close3d()
  
  # cross section
  open3d()
  plot3d(x=dat$x,y=dat$y,z=dat$z,radius=dat$radius,col=dat$color,type='s',axes=FALSE,box=FALSE,xlab='',ylab='',zlab='')
  clipplanes3d(a=0,b=1,c=0,d=0)
  umcs <- par3d()$userMatrix
  view3d(zoom=0.6,userMatrix=umcs)
  snapshot3d(filename='chrom3d_model_cross-section.png',fmt='png',width=500,height=500,webshot=FALSE)
  close3d()


  
  
  
  
  