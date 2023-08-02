
# REQUIRED PACKAGES
library(ggplot2)
library(reshape2)
library(data.table)
library(GenomicRanges)
library(maftools)
library(MutationalPatterns)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(tidyr)
library(cowplot)
library(rtracklayer)
library(RColorBrewer)



# REQUIRED FUNCTIONS
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
plot_aggregate_bigwig2 <- function(x,col='darkred',CI_color='gray80',scale_factor=1) {
  CI <- as.data.frame(t(x$CI))
  CI <- CI * scale_factor
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
cool_vioplot <- function(..., include_points=F, pt_size=1, vio_alpha=1, alpha=0.3, main=NULL, xlab=NULL, ylab=NULL, palette='Blues') {
  require(reshape2)
  require(ggplot2)
  L <- list(...)
  if(length(L)==1 & is.list(L[[1]])) {
    L <- unlist(L,recursive=F)
    nm <- names(L)
  } else {
    if(!is.null(names(L))) {
      nm <- names(L)
    } else {
      nm <- make.unique(rep('v',length(L)))
    }
  }
  names(L) <- nm
  df <- reshape2::melt(L)
  df$L1 <- factor(df$L1,levels=nm)
  p <- ggplot(df, aes(x=L1, y=value, fill=L1)) +
    geom_violin(trim=F, alpha=vio_alpha) +
    labs(title=main,x=xlab,y=ylab) 
  p <- p + scale_fill_brewer(palette=palette) + theme_classic()
  if(include_points) {
    p <- p + geom_jitter(shape=16, position=position_jitter(0.4), alpha=alpha, size=pt_size)
  }
  p + geom_boxplot(width=0.1,fill='white',outlier.shape=NA)
}
cool_sigmoid_plot2 <- function(L, special_cases=NA, stats_FUN=median, bkgd_col_low='#A1D76A', bkgd_col_high='#E9A3C9', bkgd_percent_transparency=50, pt_size=0.5) {
  require(reshape2)
  require(ggplot2)
  sigmoid = function(x) {
    1 / (1 + exp(-x))
  }
  if(!is.null(names(L))) {
    nm <- names(L)
  } else {
    nm <- make.unique(rep('v',length(L)))
  }
  names(L) <- nm
  df <- reshape2::melt(L)
  df$L1 <- factor(df$L1,levels=nm)
  if(!is.na(special_cases[1])) {
    Lu <- unname(L)
    ids <- names(unlist(Lu))
    sc <- ids %in% special_cases
    df$special_cases <- sc
  }
  X <- levels(df$L1)
  XMIN <- head(X,1)
  XMAX <- tail(X,1)
  shading=list('low'=bkgd_col_low,'high'=bkgd_col_high)
  .add_transparency <- function(color,percent=50) {
    RGB <- col2rgb(color,alpha=FALSE)
    rgb(RGB[1],RGB[2],RGB[3],max=255,alpha=(100 - percent)*255/100)
  }
  shading$low <- .add_transparency(shading$low,bkgd_percent_transparency)
  shading$high <- .add_transparency(shading$high,bkgd_percent_transparency)
  p <- ggplot(df, aes(x=L1, y=value)) +
    geom_rect(aes(xmin=XMIN,xmax=XMAX,ymin=0,ymax=Inf),fill=shading$high) +
    geom_rect(aes(xmin=XMIN,xmax=XMAX,ymin=-Inf,ymax=0),fill=shading$low)
  #uld <- unlist(dat)
  sf <- scale(sigmoid(scale(rank(unlist(L)))),scale=F,center=T) # scaling factor for nudge
  if(!is.na(special_cases[1])) {
    boolColors <- as.character(c("TRUE"="#5aae61", "FALSE"="#7b3294"))
    boolScale <- scale_colour_manual(name="special_cases", values=boolColors)
    p <- p + geom_point(aes(colour=special_cases), position=position_nudge(x=sf), size=pt_size) + boolScale + stat_summary(geom = "errorbar", fun.min = stats_FUN, fun = stats_FUN, fun.max = stats_FUN, width = .75)
  } else {
    p <- p + geom_point(position=position_nudge(x=sf), size=pt_size) + stat_summary(geom = "errorbar", fun.min = stats_FUN, fun = stats_FUN, fun.max = stats_FUN, width = .75)
  }
  p + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
            #panel.grid.major=element_blank(),
            #panel.grid.minor=element_blank(),
            panel.background=element_blank())
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
cool_plot2 <- function(x,y,bkgd='default',addTrendline=TRUE,...) {
  require('RColorBrewer')
  df <- data.frame(x,y)
  X <- densCols(x,y,colramp=colorRampPalette(c('black','white')))
  df$dens <- col2rgb(X)[1,]+1L
  cols <-  rev(colorRampPalette(brewer.pal(11,'Spectral'))(256))
  df$col <- cols[df$dens]
  #if(includeCor) {
  x[is.na(x)|is.infinite(x)] <- NA
  y[is.na(y)|is.infinite(y)] <- NA
  CorP <- paste('R=',round(cor(x,y,use='complete.obs',method='pearson'),2))
  CorS <- paste('rho=',round(cor(x,y,use='complete.obs',method='spearman'),2))
  #}
  plot(y~x,data=df[order(df$dens),],type='n',sub=paste(CorP,CorS,sep='; '),...)
  if(bkgd=='default') rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=head(cols,1)) else {
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=bkgd)
  }
  points(y~x,data=df[order(df$dens),],pch=19,col=col,...)
  #if(includeCor) {
  #  x[is.na(x)|is.infinite(x)] <- NA
  #  y[is.na(y)|is.infinite(y)] <- NA
  #  legend('topleft',paste('R=',round(cor(x,y,use='complete.obs',method=cor.method),2)),bty='n')
  #}
  if(addTrendline) {
    fit <- glm(y~x)
    co <- coef(fit)
    abline(fit, col="gray20", lwd=4, lty=2)
  }
}
cool_sigmoid_plot <- function(L, special_cases=NA, stats_FUN=median) {
  require(reshape2)
  require(ggplot2)
  sigmoid = function(x) {
    1 / (1 + exp(-x))
  }
  if(!is.null(names(L))) {
    nm <- names(L)
  } else {
    nm <- make.unique(rep('v',length(L)))
  }
  names(L) <- nm
  df <- reshape2::melt(L)
  df$L1 <- factor(df$L1,levels=nm)
  if(!is.na(special_cases[1])) {
    Lu <- unname(L)
    ids <- names(unlist(Lu))
    sc <- ids %in% special_cases
    df$special_cases <- sc
  }
  p <- ggplot(df, aes(x=L1, y=value))
  #uld <- unlist(dat)
  sf <- scale(sigmoid(scale(rank(unlist(L)))),scale=F,center=T) # scaling factor for nudge
  if(!is.na(special_cases[1])) {
    boolColors <- as.character(c("TRUE"="#5aae61", "FALSE"="#7b3294"))
    boolScale <- scale_colour_manual(name="special_cases", values=boolColors)
    p + geom_point(aes(colour=special_cases), position=position_nudge(x=sf)) + boolScale + stat_summary(geom = "errorbar", fun.min = stats_FUN, fun = stats_FUN, fun.max = stats_FUN, width = .75)
  }
  p + geom_point(position=position_nudge(x=sf)) + stat_summary(geom = "errorbar", fun.min = stats_FUN, fun = stats_FUN, fun.max = stats_FUN, width = .75)
}
sigmoid_plot <- function(L, bkgd_col_low='#A1D76A', bkgd_col_high='#E9A3C9', bkgd_alpha=0.5, pt_size=0.5, stats_FUN=median) {
  require(ggplot2)
  sigmoid = function(x) {
    1 / (1 + exp(-x))
  }
  if(!is.null(names(L))) {
    LEVELS <- names(L)
  } else {
    LEVELS <- make.unique(rep('v',length(L)))
  }
  names(L) <- LEVELS
  df <- stack(L)
  df$ind <- factor(df$ind,levels=LEVELS) # preserve original order
  df$x <- as.numeric(df$ind)
  xrng <- range(df$x)
  p <- ggplot() # initialize the plot
  # Prepare the background color
  if(!is.null(bkgd_col_low) | !is.null(bkgd_col_high)) {
    rect <- NULL
    if(!is.null(bkgd_col_low)) {
      rect <- rbind(rect, data.frame(xmin=-Inf,
                                     xmax=Inf,
                                     ymin=-Inf,
                                     ymax=0))
    }
    if(!is.null(bkgd_col_high)) {
      rect <- rbind(rect, data.frame(xmin=-Inf,
                                     xmax=Inf,
                                     ymin=0,
                                     ymax=Inf))
    }
    p <- p + geom_rect(data=rect,
                       aes(xmin=xmin,
                           xmax=xmax,
                           ymin=ymin,
                           ymax=ymax)
                       ,fill=c(bkgd_col_low,bkgd_col_high),alpha=bkgd_alpha) 
  }
  # Add the dots
  sf <- scale(sigmoid(scale(rank(unlist(L)))),scale=F,center=T) # scaling factor for nudge
  p <- p + geom_point(data=df,aes(x=x,y=values),
                      position=position_nudge(x=sf),size=pt_size)
  # Add the stats summary
  p <- p + stat_summary(data=df,mapping=aes(x=x,y=values),geom='errorbar',fun.min=stats_FUN,fun=stats_FUN,fun.max=stats_FUN)
  # Add labels/formatting
  p + scale_x_continuous(breaks=xrng[1]:xrng[2],
                         labels=LEVELS,
                         expand=c(0.01,0.01)) + 
    theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank())
}



# REQUIRED FILES
# See: https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel
# final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz (md5sum: 009b9f88ec10d07cf826174c3dea5b78)
# final_consensus_passonly.snv_mnv_indel.tcga.controlled.maf.gz (md5sum: a7f054a45e0ec4eab83cfe22b9adee73)
# simple_somatic_mutation.open.SKCA-BR.tsv.gz (md5sum: f24583ec6fb4479accde8a282ae85854)
# consensus.20170119.somatic.cna.tcga.public.tar.gz (md5sum: 0d05a8ec9f07aed1c70c6ebbe1db9aad)
# consensus.20170119.somatic.cna.icgc.public.tar.gz (md5sum: b16a44e5b52302136b2dfa2d2f4b2000)
# https://dcc.icgc.org/releases/PCAWG/consensus_sv
# 
# DSR file (this study)
icgc_maf <- '/work/devin/trad/pcawg/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz'
tcga_maf <- '/work/devin/trad/pcawg/final_consensus_passonly.snv_mnv_indel.tcga.controlled.maf.gz'
skcabr <- '/work/devin/trad/simple_somatic_mutation.open.SKCA-BR.tsv.gz' # downloaded from https://dcc.icgc.org/api/v1/download?fn=/current/Projects/SKCA-BR/simple_somatic_mutation.open.SKCA-BR.tsv.gz
icgc_cn <- '/home/devin/consensus.20170119.somatic.cna.icgc.public.tar.gz'
tcga_cn <- '/home/devin/consensus.20170119.somatic.cna.tcga.public.tar.gz'
icgc_sv <- '/home/devin/final_consensus_sv_bedpe_passonly.icgc.public.tgz'
tcga_sv <- '/home/devin/final_consensus_sv_bedpe_passonly.tcga.public.tgz'
dsr <- '/work/devin/trad/DSR_results.csv'



# Import DSR regions
chr <- c(paste0('chr',1:22),'chrX')
dsr <- read.csv(dsr)
dsr.gr <- GRanges(dsr[,1])
mcols(dsr.gr) <- dsr[,-1]
dsr.gr <- sortSeqlevels(dsr.gr)
#saveRDS(dsr.gr,'dsr.gr.RDS')



# # Import CN data from PCAWG
# untar(icgc_cn,exdir='/work/devin/trad/cn')
# untar(tcga_cn,exdir='/work/devin/trad/cn')
# cn_files <- dir('/work/devin/trad/cn',pattern='.somatic.cna.txt',full.names=TRUE)
# cn_dat <- sapply(cn_files,function(i) {
#   di <- fread(i,data.table=FALSE)
#   gri <- GRanges(seqnames=di$chromosome,
#                  ranges=IRanges(start=di$start,end=di$end),
#                  strand='*')
#   mcols(gri) <- di[,setdiff(colnames(di),c('chromosome','start','end'))]
#   seqlevelsStyle(gri) <- seqlevelsStyle(dsr.gr)
#   seqlevels(gri,pruning.mode='coarse') <- seqlevels(dsr.gr)
#   gri
# },USE.NAMES=TRUE,simplify=FALSE)
# cn_dat <- GRangesList(cn_dat)
# names(cn_dat) <- sapply(strsplit(basename(names(cn_dat)),'\\.'),'[',1)
# #saveRDS(cn_dat,'cn_dat.RDS')



# # Import PCAWG SV data
# untar(tcga_sv,exdir='/work/devin/trad/sv')
# untar(icgc_sv,exdir='/work/devin/trad/sv')
# sv_files <- dir(c('/work/devin/trad/sv/tcga/open','/work/devin/trad/sv/icgc/open'),
#                 full.names=TRUE,pattern='.somatic.sv.bedpe.gz$')
# sv_dat <- sapply(sv_files,function(i) {
#   di <- fread(i,data.table=FALSE)
#   gri <- GRanges(seqnames=di$chrom1,ranges=IRanges(start=di$start1,end=di$end1),strand=di$strand1)
#   mcols(gri) <- di[,setdiff(colnames(di),c('chrom1','start1','end1','strand1'))]
#   seqlevelsStyle(gri) <- seqlevelsStyle(dsr.gr)
#   seqlevels(gri,pruning.mode='coarse') <- seqlevels(dsr.gr)
#   gri
# },USE.NAMES=TRUE,simplify=FALSE)
# sv_dat <- GRangesList(sv_dat)
# names(sv_dat) <- sapply(strsplit(basename(names(sv_dat)),'\\.'),'[',1)
# #saveRDS(sv_dat,'sv_dat.RDS')



# Import PCAWG mutation data (ICGC and TCGA) + Brazilian cohort
COLNAMES <- c('Chromosome','Start_position','End_position','Strand','Project_Code',
              'Donor_ID','Variant_Type','Tumor_Seq_Allele1','Tumor_Seq_Allele2',
              'Tumor_Sample_Barcode')
icgc <- fread(icgc_maf)
icgc$Project_Code <- paste('ICGC',icgc$Project_Code,sep='_')
tcga <- fread(tcga_maf)
tcga$Project_Code <- paste('TCGA',tcga$Project_Code,sep='_')
skca <- icgcSimpleMutationToMAF(skcabr,MAFobj=FALSE,removeDuplicatedVariants=TRUE)
skca$Project_Code <- skca$project_code
skca$Donor_ID <- skca$icgc_donor_id
skca$Start_position <- skca$Start_Position
skca$End_position <- skca$End_Position
dat <- rbind(icgc[, ..COLNAMES],tcga[, ..COLNAMES],skca[, ..COLNAMES])



# Convert PCAWG mutations to GRanges format
gr <- GRanges(seqnames=dat$Chromosome,ranges=IRanges(start=dat$Start_position,end=dat$End_position),strand=dat$Strand,
              Project_Code=dat$Project_Code,
              Donor_ID=dat$Donor_ID,
              Variant_Type=dat$Variant_Type,
              Tumor_Seq_Allele1=dat$Tumor_Seq_Allele1,
              Tumor_Seq_Allele2=dat$Tumor_Seq_Allele2,
              Tumor_Sample_Barcode=dat$Tumor_Sample_Barcode)
seqlevelsStyle(gr) <- seqlevelsStyle(dsr.gr)
gr <- keepSeqlevels(gr,value=seqlevels(dsr.gr),pruning.mode='coarse')
gr <- sortSeqlevels(gr)
gr <- sort(gr)



# Prepare mutations object (GRanges format)
gr$mut <- paste(gr$Tumor_Seq_Allele1,gr$Tumor_Seq_Allele2,sep='>')
#saveRDS(gr,file='icgc_tcga_skca_maf.RDS')
#gr <- readRDS('/home/devin/icgc_tcga_skca_maf.RDS')

min_mutations <- 250L # required number of mutations for tumor to be included in meta-analysis

# # --- FOR DOWNSAMPLE ANALYSIS, UN-COMMENT THE FOLLOWING SECTION -------------- #
# downsample_size <- 5000L
# grls <- split(gr,gr$Tumor_Sample_Barcode)
# grls <- grls[elementNROWS(grls)>downsample_size]
# set.seed(1)
# grls <- endoapply(grls,sample,size=downsample_size,replace=FALSE)
# gr <- unlist(grls)
# min_mutations <- 100L
# # ---------------------------------------------------------------------------- #



gr.sc <- gr[gr$Project_Code %in% c('SKCA-BR','ICGC_Skin-Melanoma','TCGA_Skin-Melanoma')] # all SSM in skin cancers
gr.sc.ct <- gr.sc[gr.sc$mut %in% c('C>T','G>A')] # C>T mutations in skin cancers
gr.ct <- gr[gr$mut %in% c('C>T','G>A')] # C>T mutations in all tumors



# Mutation rate per Mb for skin tumors (DSRs)
sum(overlapsAny(gr.sc,dsr.gr[dsr.gr$DSR==TRUE])) *
  (1/length(unique(gr.sc$Tumor_Sample_Barcode))) *
  (1/sum(width(dsr.gr[dsr.gr$DSR==TRUE]))) *
  1e6 # 44.2 total mutations per Mb (in DSRs)
sum(overlapsAny(gr.sc.ct,dsr.gr[dsr.gr$DSR==TRUE])) *
  (1/length(unique(gr.sc.ct$Tumor_Sample_Barcode))) *
  (1/sum(width(dsr.gr[dsr.gr$DSR==TRUE]))) *
  1e6 # 35.9 C>T mutations per Mb (in DSRs)
sum(overlapsAny(gr.sc.ct,dsr.gr[dsr.gr$DSR==FALSE])) *
  (1/length(unique(gr.sc.ct$Tumor_Sample_Barcode))) *
  (1/sum(width(dsr.gr[dsr.gr$DSR==FALSE]))) *
  1e6 # 27.4 C>T mutations per Mb (in non-DSRs)



# For each cancer project, compute every tumor's C>T mutation rate in DSR and 
# non-DSR regions, then calculate percent difference
dsr_Mb <- sum(width(dsr.gr[dsr.gr$DSR==TRUE]))/1e6 # total bp covered by DSR regions
non_Mb <- sum(width(dsr.gr[dsr.gr$DSR!=TRUE]))/1e6 # total bp covered by non-DSR regions
# Note, regions with poor mappability have previously been removed from the DSR files
projects <- unique(gr.ct$Project_Code)
L <- sapply(projects,function(i) {
  gr.cti <- gr.ct[gr.ct$Project_Code==i] # mutations from the project
  tsb <- unique(gr.cti$Tumor_Sample_Barcode) # individual tumors
  gc()
  res <- sapply(tsb,function(ii) {
    gr.ctii <- gr.cti[gr.cti$Tumor_Sample_Barcode==ii] # mutations from the tumor
    if(length(gr.ctii) < min_mutations) NA else {
      x <- sum(overlapsAny(gr.ctii,dsr.gr[dsr.gr$DSR!=TRUE],type='within'))/non_Mb # genome-wide mutation rate per Mb (non-DSR)
      y <- sum(overlapsAny(gr.ctii,dsr.gr[dsr.gr$DSR==TRUE],type='within'))/dsr_Mb # DSR mutation rate per Mb
      c('non-DSR'=x,'DSR'=y)
    }
  },USE.NAMES=FALSE,simplify=FALSE)
  setNames(res,tsb)
},USE.NAMES=TRUE,simplify=FALSE)
#saveRDS(L,'L.RDS')
# L <- readRDS('L.RDS') # L <- readRDS('L_downsample5000.RDS') # readRDS('L_downsample.RDS')
mr <- as.data.frame(do.call(rbind,unlist(L,recursive=FALSE)))
splt <- strsplit(row.names(mr),'\\.')
mr$Project <- sapply(splt,'[',1)
mr$Tumor_Sample_Barcode <- sapply(splt,'[',2)
mr$pc <- ((mr$DSR - mr$`non-DSR`)/mr$`non-DSR`)*100 # percent-change



# Order tumors by percent-change (skin cancers only)
mrsc <- mr[mr$Project %in% c('SKCA-BR','ICGC_Skin-Melanoma','TCGA_Skin-Melanoma'),]
mrsc <- mrsc[order(mrsc$pc),]
mrsc <- mrsc[complete.cases(mrsc),]
TUMOR_ORDER <- mrsc$Tumor_Sample_Barcode
# saveRDS(mrsc,'mrsc.RDS')



# Calculate C>T fraction for each skin tumor
x <- as.data.table(mcols(gr.sc))
ct <- x[,list(sum(mut %in% c('C>T','G>A'))/.N),by=Tumor_Sample_Barcode]
ctsc <- setNames(ct$V1,ct$Tumor_Sample_Barcode)[TUMOR_ORDER]



# Calculate number of mutations in p53 pathway genes
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
p53_pathway_genes <- c('CHEK2','CHEK1','ATR','ATM','TP53','MDM4','MDM2')
p53_pathway_ids <- select(org.Hs.eg.db,
                          keys = p53_pathway_genes,
                          columns = c('SYMBOL','ENTREZID'),
                          keytype = 'SYMBOL')
gene_locs <- genes(txdb)
seqlevels(gene_locs,pruning.mode='coarse') <- seqlevels(gr)
seqinfo(gene_locs) <- seqinfo(gr)
p53_pathway_locs <- gene_locs[gene_locs$gene_id %in% p53_pathway_ids$ENTREZID]
p53_pathway_locs$gene_name <- setNames(p53_pathway_ids$SYMBOL,p53_pathway_ids$ENTREZID)[p53_pathway_locs$gene_id]
pmsc <- sapply(TUMOR_ORDER,function(i) {
  gri <- gr.sc[gr.sc$Tumor_Sample_Barcode==i]
  bai <- binnedAverage(p53_pathway_locs,coverage(gri),'x')
  bai$x <- bai$x * width(bai)
  setNames(bai$x,bai$gene_name)
},USE.NAMES=TRUE,simplify=TRUE)
pmsc <- pmsc[,TUMOR_ORDER]



# Calculate mutational signatures using MutationalPatterns R package
gr.sc$REF <- gr.sc$Tumor_Seq_Allele1
gr.sc$ALT <- gr.sc$Tumor_Seq_Allele2
gr.sc$Tumor_Sample_Barcode <- droplevels(gr.sc$Tumor_Sample_Barcode)
grl <- split(gr.sc,gr.sc$Tumor_Sample_Barcode)
snv_grl <- get_mut_type(grl,type='snv',predefined_dbs_mbs=TRUE)
genome(snv_grl) <- 'hg19'
mut_mat <- mut_matrix(vcf_list=snv_grl,ref_genome='BSgenome.Hsapiens.UCSC.hg19')
#saveRDS(mut_mat,file='mut_mat.RDS')
# mut_mat <- readRDS('mut_mat.RDS')
signatures <- get_known_signatures(muttype='snv',source='COSMIC_v3.2',genome='GRCh37')
fit_res <- fit_to_signatures(mut_mat, signatures)
sigs_to_keep <- c('SBS1','SBS2','SBS5','SBS6','SBS7a','SBS7b','SBS7c','SBS8','SBS13','SBS17a','SBS18','SBS21','SBS26','SBS32','SBS40')
strict_refit <- fit_to_signatures_strict(mut_mat, signatures[,sigs_to_keep], max_delta = 0.004)
mssc <- strict_refit$fit_res$contribution
mssc <- mssc[,TUMOR_ORDER]



# Calculate mutation rates in DSRs vs non-DSRs
DSR <- dsr.gr[dsr.gr$DSR==TRUE]
non_DSR <- dsr.gr[dsr.gr$DSR!=TRUE]
non_DSR <- resize(non_DSR,width=median(width(non_DSR)),fix='center')
# Split by TCGA, ICGC, and SKCA-BR studies
gr.sc.ct.icgc <- gr.sc.ct[gr.sc.ct$Project_Code=='ICGC_Skin-Melanoma']
gr.sc.ct.tcga <- gr.sc.ct[gr.sc.ct$Project_Code=='TCGA_Skin-Melanoma']
gr.sc.ct.skca <- gr.sc.ct[gr.sc.ct$Project_Code=='SKCA-BR']
# Mutation rate (per Mb) in DSRs
dsr.icgc <- binnedAverage(DSR,coverage(gr.sc.ct.icgc),'x')$x *
  (1/length(unique(gr.sc.ct.icgc$Tumor_Sample_Barcode))) * 1e6
non.icgc <- binnedAverage(non_DSR,coverage(gr.sc.ct.icgc),'x')$x *
  (1/length(unique(gr.sc.ct.icgc$Tumor_Sample_Barcode))) * 1e6
dsr.tcga <- binnedAverage(DSR,coverage(gr.sc.ct.tcga),'x')$x *
  (1/length(unique(gr.sc.ct.tcga$Tumor_Sample_Barcode))) * 1e6
non.tcga <- binnedAverage(non_DSR,coverage(gr.sc.ct.tcga),'x')$x *
  (1/length(unique(gr.sc.ct.tcga$Tumor_Sample_Barcode))) * 1e6
dsr.skca <- binnedAverage(DSR,coverage(gr.sc.ct.skca),'x')$x *
  (1/length(unique(gr.sc.ct.skca$Tumor_Sample_Barcode))) * 1e6
non.skca <- binnedAverage(non_DSR,coverage(gr.sc.ct.skca),'x')$x *
  (1/length(unique(gr.sc.ct.skca$Tumor_Sample_Barcode))) * 1e6



# Plot DSR v non-DSR mutation rate violin plots
options(scipen=10L)
#pdf('icgc_vio.pdf',width=2,height=2)
cool_vioplot('non-DSR'=non.icgc,'DSR'=dsr.icgc) + scale_y_log10() +
  ggtitle('ICGC') +
  ylab('C>T mutation rate\n(mutations per Mb)') + 
  theme(legend.position='none')
#dev.off()
#pdf('tcga_vio.pdf',width=2,height=2)
cool_vioplot('non-DSR'=non.tcga,'DSR'=dsr.tcga) + scale_y_log10() +
  ggtitle('TCGA') +
  ylab('C>T mutation rate\n(mutations per Mb)') + 
  theme(legend.position='none')
#dev.off()
#pdf('skca_vio.pdf',width=2,height=2)
cool_vioplot('non-DSR'=non.skca,'DSR'=dsr.skca) + scale_y_log10() +
  ggtitle('SKCA-BR') +
  ylab('C>T mutation rate\n(mutations per Mb)') + 
  theme(legend.position='none')
#dev.off()



# Plot the distribution of %-change values for skin cancers
#pdf(file='pc_vio.pdf',width=1,height=2)
cool_vioplot(mrsc$pc,include_points=TRUE,pt_size=0.2,alpha=0.8) + 
  scale_fill_brewer(palette='Greens') +
  geom_hline(yintercept=0,linetype='dashed') +
  ggtitle(label='DSR vs non-DSR\nmutation rate\nby tumor',subtitle=paste0('n=',nrow(mrsc),' tumors')) +
  labs(y='Percent change') +
  theme(legend.position='none',
        axis.text.x=element_blank())
# alternative: cool_vioplot(pc[c('SKCA-BR','ICGC_Skin-Melanoma','TCGA_Skin-Melanoma')])
#dev.off()



# Plot mutation rates around DSRs
bin_size <- 10000L
flanking <- 400e3L
profile.ctot.dsr <- aggregate_bigwig2(DSR,coverage(gr.sc.ct),bin_size=bin_size,flanking=flanking)
set.seed(1)
profile.ctot.non <- aggregate_bigwig2(sample(non_DSR,length(DSR)),coverage(gr.sc.ct),bin_size=bin_size,flanking=flanking)

profile.other.dsr <- aggregate_bigwig2(DSR,coverage(gr.sc[!(gr.sc$mut %in% c('C>T','G>A'))]),bin_size=bin_size,flanking=flanking)

profile_scale_factor <- 1/(length(unique(gr.sc.ct$Tumor_Sample_Barcode))) * (1e6/bin_size)

profile.cctt.dsr <- aggregate_bigwig2(DSR,coverage(gr.sc[gr.sc$mut %in% c('CC>TT','GG>AA')]),bin_size=bin_size,flanking=flanking)
set.seed(2)
profile.cctt.non <- aggregate_bigwig2(sample(non_DSR,length(DSR)),coverage(gr.sc[gr.sc$mut %in% c('CC>TT','GG>AA')]),bin_size=bin_size,flanking=flanking)

#pdf(file='dsr_prof.pdf',width=1.6,height=1.2)
plot_aggregate_bigwig2(profile.ctot.dsr,scale_factor=profile_scale_factor) +
  expand_limits(y=c(30,40))
#dev.off()
#pdf(file='non_prof.pdf',width=1.6,height=1.2)
plot_aggregate_bigwig2(profile.ctot.non,scale_factor=profile_scale_factor) +
  expand_limits(y=c(30,40))
#dev.off()
#pdf(file='dsr_other.pdf',width=1.6,height=1.2)
purple <- brewer.pal(3,'Purples')[3]
plot_aggregate_bigwig2(profile.other.dsr,scale_factor=profile_scale_factor,col=purple) +
  expand_limits(y=c(10,20))
#dev.off()

plot_aggregate_bigwig2(profile.cctt.dsr,scale_factor=profile_scale_factor) +
  expand_limits(y=c(0.5,0.7))
plot_aggregate_bigwig2(profile.cctt.non,scale_factor=profile_scale_factor) +
  expand_limits(y=c(0.5,0.7))



# Plot correlation between p53/ctl log2-FC and mut-rate
num_tumors <- length(unique(gr.sc.ct$Tumor_Sample_Barcode))
dsr.gr.mut <- binnedAverage(dsr.gr,coverage(gr.sc.ct),'x') # average mutations
dsr.gr.mut$mut_sum <- dsr.gr.mut$x * width(dsr.gr.mut) # sum
dsr.gr.mut$mut_rate <- dsr.gr.mut$mut_sum * (1/num_tumors) * (1/width(dsr.gr.mut)) * 1e6
# mut_rate is C>T mutations per Mb per tumor
#pdf('scatter.pdf',height=2.3333,width=2)
cool_plot2(x=dsr.gr.mut$minus_log10_pval,y=dsr.gr.mut$mut_rate,bkgd='white',addTrendline=FALSE,xlab='-log10 P-val',
           ylab='C>T mutation rate')
#dev.off()



# Dot plot, mutational signatures, and detailed tumor characterization
dat <- pivot_longer(mrsc,cols=c('non-DSR','DSR'))
dat$plot_value <- NA
non_DSR_mutation_rate <- dat[dat$name=='non-DSR',]$value
non_DSR_mutation_rate_transformed <- log10(non_DSR_mutation_rate^10)
dat[dat$name=='non-DSR',]$plot_value <- non_DSR_mutation_rate_transformed
dat[dat$name=='DSR',]$plot_value <- dat[dat$name=='DSR',]$pc + non_DSR_mutation_rate_transformed
dat$Tumor_Sample_Barcode <- factor(dat$Tumor_Sample_Barcode,levels=TUMOR_ORDER)
dat$direction <- ifelse(dat$pc >0,'pos','neg')
# Second axis
BREAKS.orig <- c(0.01,0.1,1,10,100,1000)
BREAKS <- log10(BREAKS.orig^10)
LABELS <- as.character(BREAKS.orig)
# IMPORTANT NOTE: the following dot plot contains TWO Y AXES.
# The first indicates PERCENT CHANGE in mutation rate between DSR and non-DSR
# regions for the indicated tumor (black line to pink dot). The length of the 
# black line is directly proportional to the magnitude of percent change. The
# second axis represents the BASELINE (non-DSR) mutation rate of the tumor
# (green dots). The mutation rates (green dots) undergo the transformation
# y = log10(x^10) to better visualize tumors with low overall mutation rates
# (e.g. < 1 mutation / Mb) with relatively high percent-changes 
p.dot <- ggplot(dat,aes(x=Tumor_Sample_Barcode)) +
  geom_line(aes(y=plot_value,group=Tumor_Sample_Barcode,col=direction)) +
  scale_color_manual(values=c('pos'='#404040','neg'='#2171B5')) +
  geom_point(aes(y=plot_value,fill=name),shape=21,color='white') +
  scale_fill_manual(values=c('non-DSR'='#74C476','DSR'='white')) +
  #scale_fill_manual(values=c('non-DSR'='#A1D76A','DSR'='#E9A3C9')) +
  scale_y_continuous(sec.axis=dup_axis(
    breaks=BREAKS,
    labels=LABELS)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())
# dot <- pivot_longer(mrsc,cols=c('non-DSR','DSR'))
# dot$Tumor_Sample_Barcode <- factor(dot$Tumor_Sample_Barcode,levels=TUMOR_ORDER)
# p.dot <- ggplot(dot,aes(Tumor_Sample_Barcode,value)) +
#            geom_line(aes(group=Tumor_Sample_Barcode)) +
#            geom_point(aes(color=name)) +
#            scale_y_log10() +
#            scale_color_manual(values=c('non-DSR'='#A1D76A','DSR'='#E9A3C9')) +
#            theme_bw() +
#            theme(panel.grid.major=element_blank(),
#                  panel.grid.minor=element_blank(),
#                  axis.text.x=element_blank(),
#                  axis.title.x=element_blank())

# Percent Change
pc.dat <- data.frame(x=mrsc$Tumor_Sample_Barcode,val=mrsc$pc)
pc.dat$x <- factor(pc.dat$x,levels=TUMOR_ORDER)
pc.limits <- max(abs(pc.dat$val))*(c(-1,1))
p.pc <-  ggplot(pc.dat,aes(x=x,y=1)) +
  geom_tile(aes(fill=val)) +
  #scale_fill_distiller(palette='RdGy',limits=pc.limits,direction=1) +
  scale_fill_gradient2(low='#2171B5',high='#404040',limits=pc.limits) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
# p53 pathway mutations
pm.dat <- as.data.frame.table(pmsc)
pm.dat$Var2 <- factor(pm.dat$Var2,levels=TUMOR_ORDER)
p.pm <- ggplot(pm.dat,aes(x=Var2,y=Var1,fill=Freq)) +
  geom_tile() +
  scale_fill_gradientn(colors = c('white',brewer.pal(9,'Reds')),
                       limit = c(0,6),oob=scales::squish) +
  #scale_fill_steps(limits=c(0,4),right=FALSE) +
  #scale_fill_distiller(palette='Reds',limits=c(0,4),direction=1) +
  #scale_fill_distiller(breaks=seq(-5,5,length.out=10)) +
  theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank())
# Fraction C>T mutations
ct.dat <- data.frame(x=names(ctsc),CtoT=ctsc)
ct.dat$x <- factor(ct.dat$x,levels=TUMOR_ORDER)
p.ct <- ggplot(ct.dat,aes(x=x,y=1,fill=CtoT)) +
  geom_tile() +
  scale_fill_gradientn(colors=brewer.pal(9,'Purples'),
                       limit=c(0,1)) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
# Mutational signatures
ms.dat <- as.data.frame.table(mssc)
ms.dat$Var2 <- factor(ms.dat$Var2,levels=TUMOR_ORDER)
p.ms <- ggplot(ms.dat,aes(x=Var2,y=Freq,fill=Var1)) +
  geom_bar(position='fill',stat='identity',width=1) +
  scale_fill_manual(values=c('SBS1'='#A6CEE3','SBS2'='#6A3E98','SBS5'='#C9B3D6','SBS6'='#1A79B4','SBS7a'='#FAF39D','SBS7b'='#FFCC66','SBS7c'='#F58022',
                             'SBS7d'='darkred','SBS8'='#B69064','SBS13'='#B2D88B','SBS17a'='gray','SBS18'='red',
                             'SBS21'='#A8DDB5','SBS26'='#7BCCC4','SBS32'='#D95F02','SBS40'='#984EA3')) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
# UV mechanism derived from mutational signatures
uv.dat <- apply(mssc,2,function(i) {
  uv_frac <- sum(i[c('SBS7a','SBS7b','SBS7c')])/sum(i)
  ifelse(uv_frac > 0.6,'UV','non-UV')
})
uv.dat <- data.frame(x=names(uv.dat),val=uv.dat)
uv.dat$x <- factor(uv.dat$x,levels=TUMOR_ORDER)
p.uv <- ggplot(uv.dat,aes(x=x,y=1,fill=val)) +
  geom_tile() +
  scale_fill_manual(values=c('non-UV'='gray80','UV'='#FFCC66')) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
UV_TUMORS <- unique(uv.dat[uv.dat$val=='UV',]$x)
#saveRDS(UV_TUMORS,'UV_TUMORS.RDS')

#pdf('mut_sigs.pdf',width=8,height=4)
plot_grid(
  p.dot,
  p.pc,
  p.pm,
  p.ct,
  p.ms,
  p.uv,
  align='hv',
  ncol=1,
  rel_heights=c(4,0.6,1.6,0.6,2,0.6)
)
#dev.off()



# Pan-cancer DSR mutation rate %-change plot
mr.dat <- split(mr,mr$Project)
mr.dat.pc <- lapply(mr.dat,'[[','pc')
meds <- sapply(mr.dat.pc,median,na.rm=TRUE)
mr.dat.pc <- mr.dat.pc[order(meds)]
mr.dat.pc <- mr.dat.pc[sapply(mr.dat.pc,function(x) sum(!is.na(x))>0 )]

#pdf('sigmoid_plot.pdf',width=8,height=4)
sigmoid_plot(mr.dat.pc,bkgd_alpha=0.4,pt_size=0.1)
#dev.off()



# Damage, repair, and mutations
# mutations:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76391
scatter2 <- function(x,y,nbins) {
  require(ggplot2)
  require(hexbin)
  df <- data.frame(x,y)
  ggplot(df,aes(x=x,y=y)) +
    geom_hex(bins=nbins)+
    scale_fill_distiller(palette='Spectral') +
    theme(legend.position='none',
          panel.background=element_blank())
}
dsr.gr <- readRDS('/home/devin/dsr.gr.RDS')
xr <- readRDS('/home/devin/xr.RDS')
gr <- readRDS('/home/devin/icgc_tcga_skca_maf.RDS')
#UV_TUMORS <- readRDS('/home/devin/UV_TUMORS.RDS')
#gr <- gr[gr$Tumor_Sample_Barcode %in% UV_TUMORS]

gr.sc <- gr[gr$Project_Code %in% c('TCGA_Skin-Melanoma','ICGC_Skin-Melanoma','SKCA-BR')]
gr.sc.ct <- gr.sc[gr.sc$mut %in% c('C>T','G>A')]

repair.early <- xr$CPD_1h + xr$CPD_4h + xr$CPD_8h
repair.mid <- xr$CPD_16h + xr$CPD_1d
repair.late <- xr$CPD_2d
repair.cumulative <- repair.early + repair.mid + repair.late

mutations.all <- coverage(gr.sc)
mutations.ctot <- coverage(gr.sc.ct)

dsr.gr <- binnedAverage(dsr.gr,repair.early,varname='repair.early')
dsr.gr <- binnedAverage(dsr.gr,repair.mid,varname='repair.mid')
dsr.gr <- binnedAverage(dsr.gr,repair.late,varname='repair.late')
dsr.gr <- binnedAverage(dsr.gr,repair.cumulative,varname='repair.cumulative')
dsr.gr <- binnedAverage(dsr.gr,mutations.all,varname='mutations.all')
dsr.gr <- binnedAverage(dsr.gr,mutations.ctot,varname='mutations.ctot')

dsr.gr <- dsr.gr[dsr.gr$baseMean > 0 & dsr.gr$repair.early > 0 & dsr.gr$mutations.ctot > 0 & dsr.gr$repair.late > 0]

# calculate mutation rates for all SSM and for C>T mutations
dsr.gr$mutations.all <- dsr.gr$mutations.all * width(dsr.gr) *
  (1/length(unique(gr.sc$Tumor_Sample_Barcode))) * (1e6 / width(dsr.gr))
dsr.gr$mutations.ctot <- dsr.gr$mutations.ctot * width(dsr.gr) *
  (1/length(unique(gr.sc$Tumor_Sample_Barcode))) * (1e6 / width(dsr.gr))

#pdf(file='damage.pdf',width=2,height=2)
cor.damage <- round(cor(dsr.gr$baseMean,log10(dsr.gr$mutations.ctot),method='spearman'),2)
scatter2(dsr.gr$baseMean,dsr.gr$mutations.ctot,nbins=100) +
  scale_y_log10(limits=c(1,150)) +
  xlim(c(250,1250)) +
  labs(x='Cumulative damage',
       y='C>T mutation rate\n(mut per Mb)') +
  ggtitle(
    label='Cumulative damage',
    subtitle=expr(paste(rho,'=',!!cor.damage))
  )
#dev.off()

#pdf(file='repair.pdf',width=2,height=2)
cor.repair <- round(cor(dsr.gr$repair.cumulative,log10(dsr.gr$mutations.ctot),method='spearman'),2)
scatter2(dsr.gr$repair.cumulative,dsr.gr$mutations.ctot,nbins=100) +
  scale_y_log10(limits=c(1,150)) +
  labs(x='Cumulative repair',
       y='C>T mutation rate\n(mut per Mb') +
  ggtitle(
    label='Cumulative repair',
    subtitle=expr(paste(rho,'=',!!cor.repair))
  )
#dev.off()



# # Mut p53 molecular signature (PCAWG SKCM-US)
# ### CDC20, PLK1, CENPA, KIF2C
# ### ENSG00000117399.9, ENSG00000166851.10, ENSG00000115163.10, ENSG00000142945.8
# # # REFERENCE:
# # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7546539/
# # DCC PCAWG Gene Expression
# # Gene expression file downloaded from:
# # https://dcc.icgc.org/releases/PCAWG/transcriptome/gene_expression
# # PCAWG sample sheet downloaded from: 
# # https://dcc.icgc.org/releases/PCAWG/donors_and_biospecimens
# # https://dcc.icgc.org/releases/current/Projects/SKCM-US
# TUMOR_ORDER <- readRDS('/home/devin/TUMOR_ORDER.RDS')
# mrsc <- readRDS('/home/devin/mrsc.RDS')
# ge <- fread('/home/devin/tophat_star_fpkm_uq.v2_aliquot_gl.tsv.gz',data.table=FALSE)
# pcawg_samples <- fread('/home/devin/pcawg_sample_sheet.tsv',data.table=FALSE)
# skcm_samples <- pcawg_samples[pcawg_samples$dcc_project_code=='SKCM-US',]
# skcm_rna_ids <- intersect(skcm_samples$aliquot_id,colnames(ge))
# skcm_wgs_ids <- sapply(skcm_rna_ids,function(i) {
#   spid <- skcm_samples[skcm_samples$aliquot_id==i,]$icgc_specimen_id
#   wi <- skcm_samples[skcm_samples$icgc_specimen_id==spid & skcm_samples$library_strategy=='WGS',]$aliquot_id
#   data.frame('rna'=i,'wgs'=wi)
# },USE.NAMES=FALSE,simplify=FALSE)
# skcm_wgs_ids <- do.call(rbind,skcm_wgs_ids)
# row.names(ge) <- ge$feature
# ge$feature <- NULL
# gesk <- ge[,colnames(ge) %in% skcm_rna_ids]
# colnames(gesk) <- skcm_wgs_ids[match(colnames(gesk),skcm_wgs_ids$rna),]$wgs
# # genes for mut-p53 signature
# sgenes <- c('ENSG00000117399.9',
#             'ENSG00000166851.10',
#             'ENSG00000115163.10',
#             'ENSG00000142945.8')
# sgenes <- setNames(sgenes,c(
#  'CDC20',
#  'PLK1',
#  'CENPA',
#  'KIF2C'
# ))
# 
# sgenes <- c('ENSG00000167755.9',
#             'ENSG00000261272.1',
#             'ENSG00000126545.9')
# sgenes <- setNames(sgenes,c(
#   'KLK6',
#   'MUC22',
#   'CSN1S1'))
# 
# mat <- gesk[sgenes,]
# rnk <- as.data.frame(apply(mat,1,rank))
# rnk$mut_sig <- rowMeans(rnk)
# 
# dat <- merge(x=mrsc,y=rnk,by.x='Tumor_Sample_Barcode',by.y='row.names')
# 
# 
# dat2 <- merge(x=mrsc,y=t(gesk),by.x='Tumor_Sample_Barcode',by.y='row.names')
# row.names(dat2) <- dat2$Tumor_Sample_Barcode
# dat2 <- dat2[,5:ncol(dat2)]
# corMat <- cor(dat2)





