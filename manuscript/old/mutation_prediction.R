
  setwd('~/Desktop/trad/mutation_prediction/')
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(rtracklayer)
  library(Biostrings)
  
  
  source('~/Documents/GitHub/tradtools/manuscript/old/aggregate_bigwig.R')
  
  # mutation predictions downloaded from
  # https://www.nature.com/articles/s41587-022-01353-8
  # https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-022-01353-8/MediaObjects/41587_2022_1353_MOESM4_ESM.zip
  
  kfold_data <- '~/Desktop/trad/mutation_prediction/kfold_results/Skin-Melanoma/' # https://www.nature.com/articles/s41587-022-01353-8
  blacklist <- '~/Desktop/trad/PCAWG/manuscript/MISC/hg19-blacklist.v2.bed.gz'
  mappability <- '~/Desktop/trad/PCAWG/manuscript/MISC/wgEncodeCrgMapabilityAlign100mer.bigWig'
  mutations <- '~/Desktop/trad/PCAWG/manuscript/RDS/maf_gr_03262023.RDS'
  xrseq_dir <- '~/Desktop/trad/xrseq/xrseq'
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
  xrseq_dir <- '~/Desktop/trad/xrseq/xrseq'
  gen <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  chr_to_keep <- seqlevels(gen)[1:23]
  
  
  ### Process
  # Blacklist
  bl <- read.table(blacklist,sep='\t',header=FALSE)
  bl <- GRanges(seqnames=bl[,1],ranges=IRanges(start=bl[,2]+1,end=bl[,3]))
  bl <- force_seqinfo(bl,reference=gen,chromosomes_to_keep=chr_to_keep)
  bl <- sort(bl)
  # saveRDS(bl,'bl.RDS')
  
  
  xr_ids <- c('CPD_1h','CPD_4h','CPD_8h','CPD_16h','CPD_1d','CPD_2d')
  xr_files <- dir(xrseq_dir,pattern='NHF1CPD',full.names=TRUE)
  xr <- sapply(xr_files,function(i) {
    xri <- import.bw(i)
    xri <- force_seqinfo(xri,reference=gen,chromosomes_to_keep=chr_to_keep)
    xri <- xri[!overlapsAny(xri,bl)]
    xri <- sort(xri)
    xri <- coverage(xri,weight='score')
    gc()
    xri
  },simplify=FALSE,USE.NAMES=TRUE)
  idx <- sapply(xr_ids,grep,x=names(xr),simplify=FALSE,USE.NAMES=TRUE)
  xr <- sapply(idx,function(i) {
    xi <- xr[i]
    Reduce(`+`,xi)
  },simplify=FALSE,USE.NAMES=TRUE)
  # saveRDS(xr,'xr.RDS')
  
  #xr <- readRDS('xr.RDS')
  
  # Import top DSR regions
  chr <- c(paste0('chr',1:22),'chrX')
  dsr <- read.csv('~/Desktop/trad/DSR_results.csv')
  dsr.gr <- GRanges(dsr[,1])
  mcols(dsr.gr) <- dsr[,-1]
  top_dsr <- dsr.gr[dsr.gr$DSR==TRUE]
  
  
  # Define parameters
  dsr_size <- 50e3
  bin_size <- 10e3
  flanking <- 500e3
  
  
  # Get repair rate at DSRs
  #dsr_repair <- sapply(xr,function(i) {
  #  aggregate_bigwig2(bed_file=top_dsr,bigwig_file=i,bin_size=bin_size,flanking=flanking)
  #},USE.NAMES=TRUE,simplify=FALSE)
  
  
  
  # Repair (early)
  repair_early <- xr$CPD_1h + xr$CPD_4h + xr$CPD_8h
  dsr_repair_early <- aggregate_bigwig2(top_dsr,repair_early,bin_size=bin_size,flanking=flanking)
  #pdf('dsr_repair_early.pdf',width=2,height=1)
    plot_aggregate_bigwig2(dsr_repair_early) + expand_limits(y=c(85000,115000))
  #dev.off()
  
  # Repair (late)
  repair_late <- xr$CPD_1d + xr$CPD_2d
  dsr_repair_late <- aggregate_bigwig2(top_dsr,repair_late,bin_size=bin_size,flanking=flanking)
  #pdf('dsr_repair_late.pdf',width=4,height=1)
    plot_aggregate_bigwig2(dsr_repair_late)
  #dev.off()
    
  dsr_CPD_2d <- aggregate_bigwig2(top_dsr,xr$CPD_2d,bin_size=bin_size,flanking=flanking)
  #pdf('dsr_CPD_2d.pdf',width=2,height=1)
    plot_aggregate_bigwig2(dsr_CPD_2d) + expand_limits(y=c(41000,48000))
  #dev.off()
  
    
    #pdf(file='LMNB1_profile.pdf',width=2,height=1)
    plot_aggregate_bigwig2(lamin$LMNB1,col=COLORS_DISC[5]) + expand_limits(y=c(5400,7500))
    #dev.off()
    plot_aggregate_bigwig2_heatmap(lamin$LMNB1,breaks_type='sequence',scale_factor=2)
    
  
  get_lambda_CPD_rate <- function(bed_files,fasta_file) {
    require(Biostrings)
    require(GenomicRanges)
    require(rtracklayer)
    bed <- lapply(bed_files,function(bf) {
      import.bed(bf)
    })
    bed <- do.call('c',bed)
    bed <- sort(bed)
    fa <- readDNAStringSet(fasta_file)
    dinucs <- t(dinucleotideFrequency(fa))
    colnames(dinucs) <- 'dinucs'
    cpds <- as.matrix(table(bed$name))
    colnames(cpds) <- 'cpds'
    x <- merge(x=cpds,y=dinucs,by=0,)
    setNames(x$cpds / x$dinucs,nm=x$Row.names)
  }
  
  plot2 <- function(x,y,xlab=deparse(substitute(x)),ylab=deparse(substitute(y)),log_transf=TRUE,cor_method='pearson',...) {
    require(ggplot2)
    require(ggpointdensity)
    stopifnot(length(x)==length(y))
    stat <- cor(x,y,method=cor_method)
    stat <- paste('cor =',signif(stat,3))
    data <- data.frame(x,y)
    p <- ggplot(data,aes(x=x,y=y)) +
      #stat_density_2d(aes(fill=after_stat(density)),geom='raster',contour=FALSE) +
      #scale_fill_distiller(palette='Spectral',direction=-1) +
      geom_pointdensity(...) +
      scale_color_distiller(palette='Spectral') +
      xlab(xlab) +
      ylab(ylab) +
      theme_bw() +
      theme(legend.position='none') +
      annotate(geom='text',label=stat,x=Inf,y=Inf,hjust=1,vjust=1)
    if(log_transf) {
      p <- p + scale_x_log10(expand=c(0,0)) + scale_y_log10(expand=c(0,0))
    }
    #p + geom_smooth()
    p
  }
  
  mappable_sequence_counts <- function(gr,genome,mappability_track,sequences,chunk_size,prefix='mappable_',verbose=T) {
    t0 <- Sys.time()
    require(GenomicRanges)
    require(data.table)
    if(verbose) message('Counting each occurrence of the following sequences in granges: ',paste(sequences,collapse=', '))
    stopifnot(is(gr,'GenomicRanges'))
    stopifnot(is(genome,'BSgenome'))
    if(!(all(seqlevels(gr) %in% names(genome)))) stop('One or more seqlevels in "gr" not found in "genome"')
    stopifnot(is(mappability_track,'GenomicRanges'))
    stopifnot(is(sequences,'character'))
    stopifnot(isSingleNumber(chunk_size))
    num_chars <- nchar(sequences)
    stopifnot(all(num_chars == num_chars[1]))
    dat <- granges(gr)
    wi <- num_chars[1]
    nr <- length(dat)
    nc <- length(sequences)
    start_idx <- seq(1,length(dat),by=chunk_size)
    end_idx <- start_idx + chunk_size - 1
    end_idx[length(end_idx)] <- length(dat)
    num_chunks <- length(start_idx)
    nuc_counts <- data.table(matrix(data=NA_integer_,nrow=nr,ncol=nc))
    colnames(nuc_counts) <- sequences
    new_colnames <- paste0(prefix,colnames(nuc_counts))
    for(i in seq_along(start_idx)) {
      if(verbose) message('\rProcessing chunk ',i,' of ',num_chunks,appendLF=i == num_chunks)
      idx <- start_idx[i]:end_idx[i]
      dati <- dat[idx]
      ID <- seq_along(dati)
      resi <- data.table(ID=ID)
      dati$ID <- ID
      opi <- findOverlapPairs(dati,mappability_track)
      pii <- pintersect(opi)
      bsgv <- Views(genome,pii)
      cti <- oligonucleotideFrequency(x=bsgv,width=wi,with.labels=T)
      cti <- data.table(cbind(cti,ID=S4Vectors::first(opi)$ID))
      cti <- cti[,lapply(.SD,sum),by=ID,.SDcols=sequences]
      resi <- data.table::merge.data.table(resi,cti,by='ID',all.x=T)
      nuc_counts[idx] <- resi[,-1]
    }
    colnames(nuc_counts) <- new_colnames
    if(any(idx_to_remove <- new_colnames %in% colnames(mcols(gr)))) {
      if(verbose) message('column(s) ',paste(new_colnames[idx_to_remove],collapse=', '),' already exist and will be overwritten')
      mcols(gr)[new_colnames[idx_to_remove]] <- NULL
    }
    mcols(gr) <- cbind(mcols(gr),nuc_counts)
    tf <- Sys.time()-t0
    if(verbose) message('Done. Elapsed time: ',format(round(tf,2)))
    return(gr)    
  }
  
  mappable_total <- function(gr,mappability_track,verbose=TRUE) {
    require('GenomicRanges')
    stopifnot(is(gr,'GRanges'))
    stopifnot(is(mappability_track,'GRanges'))
    stopifnot(isTRUEorFALSE(verbose))
    if(any(colnames(mcols(gr))=='mappable_total')) {
      if(verbose) message('column "mappable_total" already exists and will be overwritten')
      mcols(gr)['mappable_total'] <- NULL
    }
    if(verbose) message('Calculating total number of mappable bases in granges ...',appendLF=F)
    mtc <- coverage(mappability_track)
    mtc <- mtc[seqlevels(gr)]
    gr$mappable_total <- GenomicRanges::binnedAverage(bins=gr,numvar=mtc,'x')$x * width(gr)
    message(' DONE.')
    return(gr)
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
  
  binnedSum <- function (bins, numvar, varname, na.rm = FALSE) 
  {
    if (!is(bins, "GRanges")) 
      stop("'x' must be a GRanges object")
    if (!is(numvar, "RleList")) 
      stop("'numvar' must be an RleList object")
    if (!identical(seqlevels(bins), names(numvar))) 
      stop("'seqlevels(bin)' and 'names(numvar)' must be identical")
    viewMeans2 <- function(v, na.rm = FALSE) {
      if (!isTRUEorFALSE(na.rm)) 
        stop("'na.rm' must be TRUE or FALSE")
      means <- viewMeans(v, na.rm = na.rm)
      w0 <- width(v)
      v1 <- trim(v)
      w1 <- width(v1)
      if (na.rm) {
        na_count <- sum(is.na(v1))
        w0 <- w0 - na_count
        w1 <- w1 - na_count
      }
      means <- means * w1/w0
      means[w0 != 0L & w1 == 0L] <- 0
      means
    }
    bins_per_chrom <- split(ranges(bins), seqnames(bins))
    means_list <- lapply(names(numvar), function(seqname) {
      v <- Views(numvar[[seqname]], bins_per_chrom[[seqname]])
      viewMeans2(v, na.rm = na.rm)
    })
    new_mcol <- unsplit(means_list, as.factor(seqnames(bins))) * width(bins)
    mcols(bins)[[varname]] <- new_mcol
    bins
  }
  
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
  
  import_kfold <- function(data_dir,genome,chromosomes_to_keep) {
    require('GenomicRanges')
    fls <- dir(path=data_dir,pattern='.txt',full.names=TRUE)
    grs <- lapply(fls,function(i) {
      src <- basename(tools::file_path_sans_ext(i))
      gri <- read.table(i,sep='\t',header=TRUE)
      mc <- cbind(gri[,-c(1,2,3)],'source'=src)
      gr <- GRanges(seqnames=gri[,'CHROM'],ranges=IRanges(start=gri[,'START']+1,end=gri[,'END']))
      mcols(gr) <- mc
      force_seqinfo(gr,reference=genome,chromosomes_to_keep=chromosomes_to_keep)
    })
    gr <- do.call('c',grs)
    gr <- sort(gr)
  }
  
  
  # Import data
  kfold_data <- '~/Desktop/trad/mutation_prediction/kfold_results/Skin-Melanoma/' # https://www.nature.com/articles/s41587-022-01353-8
  blacklist <- '~/Desktop/trad/PCAWG/manuscript/MISC/hg19-blacklist.v2.bed.gz'
  mappability <- '~/Desktop/trad/PCAWG/manuscript/MISC/wgEncodeCrgMapabilityAlign100mer.bigWig'
  mutations <- '~/Desktop/trad/PCAWG/manuscript/RDS/maf_gr_03262023.RDS'
  xrseq_dir <- '~/Desktop/trad/xrseq/xrseq'
  #cpds.1 <- '~/Desktop/trad/results/bw/IMR_200J_UVB_1_possort_markdup_spec1.bw'
  #cpds.2 <- '~/Desktop/trad/results/bw/IMR_200J_UVB_2_possort_markdup_spec1.bw'
  
  #cpds.1 <- '~/Desktop/trad/results/bw/IMR_100J_UVC_1_possort_markdup_spec1.bw'
  #cpds.2 <- '~/Desktop/trad/results/bw/IMR_100J_UVC_2_possort_markdup_spec1.bw'
  
  gen <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  chr_to_keep <- seqlevels(gen)[1:22]
  
  ### Process
  # Blacklist
  bl <- read.table(blacklist,sep='\t',header=FALSE)
  bl <- GRanges(seqnames=bl[,1],ranges=IRanges(start=bl[,2]+1,end=bl[,3]))
  bl <- force_seqinfo(bl,reference=gen,chromosomes_to_keep=chr_to_keep)
  bl <- sort(bl)
  # saveRDS(bl,'bl.RDS')
  
  # Mappability
  mt <- import.bw(mappability)
  mt <- force_seqinfo(mt,reference=gen,chromosomes_to_keep=chr_to_keep)
  mt <- sort(mt)
  # saveRDS(mt,'mt.RDS')
  
  # Selection
  pred <- import_kfold(kfold_data,genome=gen,chromosomes_to_keep=chr_to_keep) 
  pred <- pred[pred$source != 'Regions_withheld_from_all_folds']
  anyDuplicated(pred)
  pred <- pred[!overlapsAny(pred,bl)]
  # saveRDS(pred,'pred.RDS')
  
  # C>T melanoma mutations, expanded cohort
  mut <- readRDS(mutations)
  mut <- force_seqinfo(mut,reference=gen,chromosomes_to_keep=chr_to_keep)
  mut <- mut[!overlapsAny(mut,bl)] # remove mutations overlapping blacklisted loci
  mutc <- coverage(mut)
  # saveRDS(mutc,'mutc.RDS')
  
  # Repair
  xr_ids <- c('CPD_1h','CPD_4h','CPD_8h','CPD_16h','CPD_1d','CPD_2d')
  xr_files <- dir(xrseq_dir,pattern='NHF1CPD',full.names=TRUE)
  xr <- sapply(xr_files,function(i) {
    xri <- import.bw(i)
    xri <- force_seqinfo(xri,reference=gen,chromosomes_to_keep=chr_to_keep)
    xri <- xri[!overlapsAny(xri,bl)]
    xri <- sort(xri)
    xri <- coverage(xri,weight='score')
    gc()
    xri
  },simplify=FALSE,USE.NAMES=TRUE)
  idx <- sapply(xr_ids,grep,x=names(xr),simplify=FALSE,USE.NAMES=TRUE)
  xr <- sapply(idx,function(i) {
    xi <- xr[i]
    Reduce(`+`,xi)
  },simplify=FALSE,USE.NAMES=TRUE)
  # saveRDS(xr,'xr.RDS')
  
  # Damage
  uvc <- list(rep1='~/Desktop/trad/results/bed/IMR_100J_UVC_1_possort_markdup_spec1.bed.gz',rep2='~/Desktop/trad/results/bed/IMR_100J_UVC_2_possort_markdup_spec1.bed.gz')
  uvb <- list(rep1='~/Desktop/trad/results/bed/IMR_200J_UVB_1_possort_markdup_spec1.bed.gz',rep2='~/Desktop/trad/results/bed/IMR_200J_UVB_2_possort_markdup_spec1.bed.gz')  
  p53 <- list(rep1='~/Desktop/trad/results/bed/TP53_200J_UVB_1_possort_markdup_spec1.bed.gz',rep2='~/Desktop/trad/results/bed/TP53_200J_UVB_2_possort_markdup_spec1.bed.gz')
  mel <- list(rep1='~/Desktop/trad/results/bed/Melanocytes_200J_UVB_1_possort_markdup_spec1.bed.gz',rep2='~/Desktop/trad/results/bed/Melanocytes_200J_UVB_2_possort_markdup_spec1.bed.gz')
  cpd_samples <- setNames(list(uvc,uvb,p53,mel),c('uvc','uvb','p53','mel'))
  cpds <- lapply(cpd_samples,function(i) {
    r1 <- import.bed(i$rep1)
    r2 <- import.bed(i$rep2)
    x <- c(r1,r2)
    x <- force_seqinfo(x,reference=gen,chromosomes_to_keep=chr_to_keep)
    x <- x[!overlapsAny(x,bl)]
    sort(x)
  })
  names(cpds) <- names(cpd_samples)
  # saveRDS(cpds,'cpds.RDS')
  cpds <- lapply(cpds,function(i) {
    list(TT=coverage(i[i$name=='TT']),
         TC=coverage(i[i$name=='TC']),
         CT=coverage(i[i$name=='CT']),
         CC=coverage(i[i$name=='CC']))
  })
  names(cpds) <- names(cpd_samples)
  # saveRDS(cpds,'cpds2.RDS')
  
  
  
  # Putting it all together
  mt <- readRDS('mt.RDS')
  pred <- readRDS('pred.RDS')
  mutc <- readRDS('mutc.RDS')
  cpds <- readRDS('cpds2.RDS')
  xr <- readRDS('xr.RDS')
  
  #dat <- pred
  
  sl <- seqlengths(gen)[chr_to_keep]
  bins <- tileGenome(sl,tilewidth=1e6,cut.last.tile.in.chrom=TRUE)
  dat <- bins
  
  
  dat <- binnedSum(dat,mutc,'CtoT')
  ul <- unlist(cpds,recursive=FALSE)
  cpd_counts <- as.data.frame(sapply(ul,function(i) {
    binnedSum(dat,i,'x')$x
  },simplify=FALSE,USE.NAMES=TRUE))
  mcols(dat) <- cbind(mcols(dat),cpd_counts)
  repair <- as.data.frame(sapply(xr,function(i) {
    binnedSum(dat,i,'x')$x
  },simplify=FALSE,USE.NAMES=TRUE))
  mcols(dat) <- cbind(mcols(dat),repair)
  dat <- mappable_sequence_counts(dat,genome=gen,mappability_track=mt,sequences=c('TT','TC','CT','CC'),chunk_size=20e3)
  dat <- mappable_sequence_counts(dat,genome=gen,mappability_track=mt,sequences=c('C','G'),chunk_size=20e3)
  dat <- binnedAverage(dat,coverage(mt,weight='score'),'mappability')
  
  
  # saveRDS(dat,'dat.RDS')
  
  
  
  
  
  
  
  
  #dat <- readRDS('dat.RDS')
  # dat <- readRDS('dat_100kb.RDS')
  dat <- readRDS('dat_1Mb.RDS')
  dat <- dat[complete.cases(mcols(dat))]
  dat <- dat[dat$mappability > 0.8]
  
  # # Pre-filter
  # num.idx <- sapply(mcols(dat),class) %in% c('numeric','integer')
  # mat <- as.matrix(mcols(dat)[,num.idx])
  # #mat <- mat[complete.cases(mat),]
  # map.idx <- grep('mappable_',colnames(mat))
  # rep.idx <- grep('CPD_',colnames(mat))
  # tra.idx <- grep('\\.[T|C][T|C]',colnames(mat))
  # map.rs <- rowSums(mat[,map.idx])
  # rep.rs <- rowSums(mat[,rep.idx])
  # tra.rs <- rowSums(mat[,tra.idx])
  # percentile <- 0.01
  # to_keep <- map.rs > quantile(map.rs,percentile,na.rm=TRUE) & rep.rs > quantile(rep.rs,percentile,na.rm=TRUE) & tra.rs > quantile(tra.rs,percentile,na.rm=TRUE)
  # dat <- dat[to_keep]
  
  

  # Positive selection
  dat$pos_selection <- dat$Y_PRED - dat$Y_TRUE
  # Mutation rate
  dat$mut_rate <- (dat$CtoT / (dat$mappable_C + dat$mappable_G))*1000
  # Susceptibility
  trad <- mcols(dat)[,grep('\\.[T|C][T|C]',colnames(mcols(dat)))]
  susc <- sapply(colnames(trad),function(cn) {
    x <- trad[,cn,drop=TRUE]
    if(grepl('\\.TT$',cn)) den <- dat$mappable_TT else {
      if(grepl('\\.TC$',cn)) den <-  dat$mappable_TC else {
        if(grepl('\\.CT$',cn)) den <- dat$mappable_CT else {
          if(grepl('\\.CC$',cn)) den <- dat$mappable_CC else {
            den <- NA
          }
        }
      }
    }
    x / den
  },simplify=TRUE,USE.NAMES=TRUE)
  colnames(susc) <- paste(colnames(susc),'seq_norm',sep='.')
  mcols(dat) <- cbind(mcols(dat),susc)
  
  # Total CPDs
  dat$TRAD_count.raw <- rowSums(as.matrix(trad))
  dat$TRAD_count.seq_norm <- dat$TRAD_count.raw / (dat$mappable_TT + dat$mappable_TC + dat$mappable_CT + dat$mappable_CC)
  dat$TRAD_count.C_raw <- rowSums(as.matrix(trad[,grep('TC|CT|CC',colnames(trad))]))
  dat$TRAD_count.C_seq_norm <- dat$TRAD_count.C_raw / (dat$mappable_TC + dat$mappable_CT + dat$mappable_CC)
  

  tmpfun <- function(x,n=10) colSums(matrix(x,nrow=n,byrow=TRUE))
  
  
  # Mutations are distributed unevenly across the cancer genome and mutation 
  # rates across genomic regions are highly heterogeneous [12] due to genomic 
  # and epigenetic features including cytosine methylation [13], replication 
  # timing [14], tri-nucleotide/penta-nucleotide context composition [5], 
  # transcription factor binding, chromatin organization [15], gene expression 
  # levels [16], orientation of the DNA minor groove around nucleosomes [17], 
  # CTCF binding [18] and gene body features such as introns and exons [19].
  
  # Susceptibility: (O - E) / E 
  # or (O - E)^2 / E
  
  # expected: 
  # uvc
  fasta_file <- '~/Desktop/tp53/lambda.fa.gz'
  
  
  uvc_expRate <- get_lambda_CPD_rate(bed_files=c('~/Desktop/trad/results/bed/IMR_100J_UVC_1_possort_markdup_spec2.bed.gz',
                                         '~/Desktop/trad/results/bed/IMR_100J_UVC_2_possort_markdup_spec2.bed.gz'),
                             fasta_file=fasta_file)
  dat$uvc_obs <- dat$uvc.TT + dat$uvc.TC + dat$uvc.CT + dat$uvc.CC
  dat$uvc_exp <- dat$mappable_TT*uvc_expRate[['TT']] +
                 dat$mappable_TC*uvc_expRate[['TC']] +
                 dat$mappable_CT*uvc_expRate[['CT']] +
                 dat$mappable_CC*uvc_expRate[['CC']]
  dat$uvc_exp <- dat$uvc_exp * (sum(dat$uvc_obs) / sum(dat$uvc_exp))
  dat$uvc_susc <- (dat$uvc_obs - dat$uvc_exp) / dat$uvc_exp
  
  
  uvb_expRate <- get_lambda_CPD_rate(bed_files=c('~/Desktop/trad/results/bed/IMR_200J_UVB_1_possort_markdup_spec2.bed.gz',
                                                 '~/Desktop/trad/results/bed/IMR_200J_UVB_2_possort_markdup_spec2.bed.gz'),
                                     fasta_file=fasta_file)
  dat$uvb_obs <- dat$uvb.TT + dat$uvb.TC + dat$uvb.CT + dat$uvb.CC
  dat$uvb_exp <- dat$mappable_TT*uvb_expRate[['TT']] +
    dat$mappable_TC*uvb_expRate[['TC']] +
    dat$mappable_CT*uvb_expRate[['CT']] +
    dat$mappable_CC*uvb_expRate[['CC']]
  dat$uvb_exp <- dat$uvb_exp * (sum(dat$uvb_obs) / sum(dat$uvb_exp))
  dat$uvb_susc <- (dat$uvb_obs - dat$uvb_exp) / dat$uvb_exp
  
  
  mel_expRate <- get_lambda_CPD_rate(bed_files=c('~/Desktop/trad/results/bed/Melanocytes_200J_UVB_1_possort_markdup_spec2.bed.gz',
                                                 '~/Desktop/trad/results/bed/Melanocytes_200J_UVB_2_possort_markdup_spec2.bed.gz'),
                                     fasta_file=fasta_file)
  dat$mel_obs <- dat$mel.TT + dat$mel.TC + dat$mel.CT + dat$mel.CC
  dat$mel_exp <- dat$mappable_TT*mel_expRate[['TT']] +
    dat$mappable_TC*mel_expRate[['TC']] +
    dat$mappable_CT*mel_expRate[['CT']] +
    dat$mappable_CC*mel_expRate[['CC']]
  dat$mel_exp <- dat$mel_exp * (sum(dat$mel_obs) / sum(dat$mel_exp))
  dat$mel_susc <- (dat$mel_obs - dat$mel_exp) / dat$mel_exp
  
  
  p53_expRate <- get_lambda_CPD_rate(bed_files=c('~/Desktop/trad/results/bed/TP53_200J_UVB_1_possort_markdup_spec2.bed.gz',
                                                 '~/Desktop/trad/results/bed/TP53_200J_UVB_2_possort_markdup_spec2.bed.gz'),
                                     fasta_file=fasta_file)
  dat$p53_obs <- dat$p53.TT + dat$p53.TC + dat$p53.CT + dat$p53.CC
  dat$p53_exp <- dat$mappable_TT*p53_expRate[['TT']] +
    dat$mappable_TC*p53_expRate[['TC']] +
    dat$mappable_CT*p53_expRate[['CT']] +
    dat$mappable_CC*p53_expRate[['CC']]
  dat$p53_exp <- dat$p53_exp * (sum(dat$p53_obs) / sum(dat$p53_exp))
  dat$p53_susc <- (dat$p53_obs - dat$p53_exp) / dat$p53_exp
  
  
  

  
  
  
  
  

  
  # Neural network 
  # https://datascience.stackexchange.com/questions/44644/how-to-determine-feature-importance-in-a-neural-network
  # https://www.analyticsvidhya.com/blog/2017/09/creating-visualizing-neural-network-in-r/
  # https://www.r-bloggers.com/2021/04/deep-neural-network-in-r/
  # https://www.learnbymarketing.com/tutorials/neural-networks-in-r-tutorial/
  
  # Study this: https://github.com/h2oai/h2o-tutorials/tree/master/tutorials/deeplearning
  # Study this: https://towardsdatascience.com/a-guide-to-using-h2o-ai-in-r-99cf6265bc05
  
  #library(neuralnet)
  library(h2o)
  # nn_dat <- cbind(CtoT=all_dat$CtoT,
  #                 pos_selection=all_dat$pos_selection,
  #                 all_dat[,grepl('mappable_',colnames(all_dat))],
  #                 all_dat[,grepl('\\.[T|C][T|C]$',colnames(all_dat))],
  #                 all_dat[,grepl('^CPD_',colnames(all_dat))])
  # 
  # nn_dat <- apply(nn_dat,2,tmpfun)
  all_dat <- mcols(dat)
  nn_dat <- cbind(mut_rate=dat$mut_rate,
                  all_dat[,grepl('\\.seq_norm$',colnames(all_dat))],
                  all_dat[,grepl('^CPD_',colnames(all_dat))])
  nn_dat <- as.data.frame(nn_dat)
  # scale data
  m <- colMeans(nn_dat)
  s <- apply(nn_dat, 2, sd)
  nn_dat <- scale(nn_dat, center = m, scale = s)
  
  # train set and test set
  idx <- sample(c(TRUE,FALSE),size=nrow(nn_dat),replace=TRUE,prob=c(0.8,0.2))
  train <- nn_dat[idx,]
  test <- nn_dat[!idx,]
  localh2o <- h2o.init(ip='localhost',port=54321,startH2O=TRUE,nthreads=6)
  train <- as.h2o(train)
  test <- as.h2o(test)
  response <- 'mut_rate'
  predictors <- setdiff(colnames(nn_dat),y)
  
  # all model comparison
  aml <- h2o.automl(x=predictors,
                    y=response,
                    training_frame=train,
                    validation_frame=test,
                    max_models=15,
                    seed=1)
  
  
  
  # prepare model
  
  
  

  
  
  model <- h2o.deeplearning(x=x,
                            y=y,
                            training_frame=train,
                            validation_frame=test,
                            distribution='AUTO',
                            hidden=c(200,200),
                            input_dropout_ratio=0.2,
                            l1=1e-5,
                            epochs=1000)
  
  model.gbm <- h2o.gbm(x=x,
                       y=y,
                       training_frame=train,
                       validation_frame=test)
  
  model.rf <- h2o.randomForest(x=x,
                               y=y,
                               training_frame=train,
                               validation_frame=test)
  
  model.glm <- h2o.glm(x=x,
                       y=y,
                       training_frame=train,
                       validation_frame=test)
  
  
  
  # n <- neuralnet(mut_rate ~ .,data=train,
  #                algorithm='rprop+',
  #                hidden=c(10,3),
  #                threshold=0.1,
  #                lifesign='full')
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  plot2(dat$CPD_1h,dat$mut_rate,cor_method='spearman')
  plot2(dat$TRAD_count.seq_norm,dat$mut_rate,cor_method='spearman')
  
  plot2(log2(dat$p53.TT + dat$p53.TC + dat$p53.CT + dat$p53.CC + 1 / dat$uvb.TT + dat$uvb.TC + dat$uvb.CT + dat$uvb.CC + 1), dat$mut_rate,cor_method='spearman')
  
  
  plot2(log2(dat$uvc.TT + dat$uvc.TC + dat$uvc.CT + dat$uvc.CC + 1 / dat$uvb.TT + dat$uvb.TC + dat$uvb.CT + dat$uvb.CC + 1), dat$mut_rate,cor_method='spearman',xlab='UV-C diff. susc',ylab='C>T mutation rate')
  
  
  plot2(dat$p53.TT + dat$p53.TC + dat$p53.CT + dat$p5)
  
  # # Intergenic mutations
  # library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  # txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  # # genic <- genes(txdb)
  # # genic <- force_seqinfo(genic,reference=gen,chromosomes_to_keep=chr_to_keep)
  # # strand(genic) <- '*'
  # # genic <- reduce(genic)
  # coding <- cds(txdb)
  # coding <- force_seqinfo(coding,reference=gen,chromosomes_to_keep=chr_to_keep)
  # strand(coding) <- '*'
  # coding <- reduce(coding)
  # mutc <- readRDS('mutc.RDS')
  # mutcig <- as(mutc,'GRanges')
  # mutcig <- mutcig[!overlapsAny(mutcig,coding)]
  # mutcig <- coverage(mutcig,weight='score')
  # 
  # 
  # dat <- binnedSum(dat,mutcig,'noncoding_CtoT')
  #dat2 <- dat[!overlapsAny(dat,coding)]
  
  
  all_dat <- as.data.frame(mcols(dat)[,sapply(mcols(dat),class) %in% c('numeric','integer')])
  #all_dat <- scale(all_dat)
  all_dat <- cbind(Y_PRED=all_dat$Y_PRED,Y_TRUE=all_dat$Y_TRUE,CtoT=all_dat$CtoT,mut_rate=all_dat$mut_rate,
                   all_dat[,grepl('\\.[T|C][T|C]\\.seq_norm',colnames(all_dat))],
                   all_dat[,grepl('^CPD_',colnames(all_dat))])
  
  cm <- cor(all_dat,method='spearman')
  
  pheatmap(cm)
  
  
  tmp <- cm[,'mut_rate']
  tmp[order(tmp,decreasing=TRUE)]
  
  
  
  
  
  # 
  # 
  # 
  # # TT_rate
  # dat$TT_rate <- dat$uvc.TT / dat$mappable_TT
  # # TC_rate
  # dat$TC_rate <- dat$uvc.TC / dat$mappable_TC
  # # CT_rate
  # dat$CT_rate <- dat$uvc.CT / dat$mappable_CT
  # # CC_rate
  # dat$CC_rate <- dat$uvc.CC / dat$mappable_CC
  # # C-dipyrimidine rate
  # dat$Cdipyr_rate <- (dat$uvc.TC + dat$uvc.CT + dat$uvc.CC) / (dat$mappable_TC + dat$mappable_CT + dat$mappable_CC)
  # 
  # 
  # mat <- mcols(dat)[,-4]
  # cm <- cor(as.matrix(mat),method='spearman')
  # 
  # 
  # binnedSum(dat,cpds$uvc$TT,'uvc_TT')
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # # # Note blacklisted counts already removed during TRADtools processing for .bw files
  # # uvc1 <- import.bw('~/Desktop/trad/results/bw/IMR_100J_UVC_1_possort_markdup_spec1.bw')
  # # uvc2 <- import.bw('~/Desktop/trad/results/bw/IMR_100J_UVC_2_possort_markdup_spec1.bw')
  # # uvc1 <- force_seqinfo(uvc1,reference=gen,chromosomes_to_keep=chr_to_keep)
  # # uvc2 <- force_seqinfo(uvc2,reference=gen,chromosomes_to_keep=chr_to_keep)
  # # uvc <- coverage(uvc1,weight='score') + coverage(uvc2,weight='score')
  # # 
  # # uvb1 <- import.bw('~/Desktop/trad/results/bw/IMR_200J_UVB_1_possort_markdup_spec1.bw')
  # # uvb2 <- import.bw('~/Desktop/trad/results/bw/IMR_200J_UVB_2_possort_markdup_spec1.bw')
  # # uvb1 <- force_seqinfo(uvb1,reference=gen,chromosomes_to_keep=chr_to_keep)
  # # uvb2 <- force_seqinfo(uvb2,reference=gen,chromosomes_to_keep=chr_to_keep)
  # # uvb <- coverage(uvb1,weight='score') + coverage(uvb2,weight='score')
  # # 
  # # tp531 <- import.bw('~/Desktop/trad/results/bw/TP53_200J_UVB_1_possort_markdup_spec1.bw')
  # # tp532 <- import.bw('~/Desktop/trad/results/bw/TP53_200J_UVB_2_possort_markdup_spec1.bw')
  # # tp531 <- force_seqinfo(tp531,reference=gen,chromosomes_to_keep=chr_to_keep)
  # # tp532 <- force_seqinfo(tp532,reference=gen,chromosomes_to_keep=chr_to_keep)
  # # tp53 <- coverage(tp531,weight='score') + coverage(tp532,weight='score')
  # # 
  # 
  # tmp <- mappable_sequence_counts(tmp,genome=gen,mappability_track=mt,sequences=c('TT','TC','CT','CC'),chunk_size=20e3)
  # 
  # 
  # tmp$susc <- (tmp$uvc) / (tmp$mappable_TT + tmp$mappable_TC + tmp$mappable_CT + tmp$mappable_CC)
  # 
  # # Repair
  # 
  # 
  # 
  # tmp <- binnedSum(pred,cpdc,varname='CPDs')
  # tmp <- binnedSum(tmp,mutc,varname='muts')
  # tmp <- binnedSum(tmp,uvc,varname='uvc')
  
