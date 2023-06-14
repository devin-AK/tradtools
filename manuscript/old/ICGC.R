
  setwd('~/Desktop/trad/PCAWG/manuscript/')
 # setwd('/work/data/trad')
  
  
  
  ### Load libraries
  library(R.utils)
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(pheatmap)
  library(RColorBrewer)
  library(maftools)
  library(sigminer)
  library(parallel)
  library(ggplot2)
  library(tidyr)
  
  
  
  # Load functions
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
  
  mutation_signatures_from_maf <- function(maf, output_file=NULL, ref_genome='BSgenome.Hsapiens.UCSC.hg19', signature_index='ALL', signature_db='SBS_hg19') {
    require('maftools')
    require('sigminer')
    if(is.character(maf)) {
      maf <- maftools::read.maf(maf,useAll=TRUE,removeDuplicatedVariants=TRUE)
    }
    stopifnot(is(maf,'MAF'))
    if(!is.null(output_file)) {
      ext <- tools::file_ext(output_file)
      if(nchar(ext)==0) output_file <- paste0(output_file,'.csv')
      stopifnot(tools::file_ext(output_file)=='csv')
    }
    # process
    st <- sigminer::sig_tally(maf,ref_genome=ref_genome,use_syn=TRUE)
    mat <- t(st$nmf_matrix)
    sf <- sigminer::sig_fit(mat,sig_index=signature_index,sig_db=signature_db)
    if(!is.null(output_file)) {
      write.csv(t(sf),file=output_file)
    } else {
      return(sf)
    }
  }
  
  convert_maf_to_granges <- function(maf, output_file=NULL) {
    require(maftools)
    require(GenomicRanges)
    require(data.table)
    if(is.character(maf)) {
      maf <- maftools::read.maf(maf,useAll=TRUE,removeDuplicatedVariants=TRUE)
    }
    stopifnot(is(maf,'MAF'))
    if(!is.null(output_file)) {
      ext <- tools::file_ext(output_file)
      if(nchar(ext)==0) output_file <- paste0(output_file,'.RDS')
      stopifnot(toupper(tools::file_ext(output_file))=='RDS')
    }
    stopifnot(identical(colnames(maf@data),colnames(maf@maf.silent)))
    maf <- rbind(maf@data,maf@maf.silent)
    stopifnot(is(maf,'data.table'))
    stopifnot(identical(maf$Reference_Allele,maf$Tumor_Seq_Allele1))
    cn <- colnames(maf)
    cols <- c('project_code','donor_id','tumor_sample_barcode','reference_allele','tumor_seq_allele2','variant_type')
    stopifnot(all(sapply(cols,function(i)sum(grepl(i,cn,ignore.case=TRUE))==1)))
    cols_to_keep <- sapply(cols,function(i)grep(i,cn,ignore.case=TRUE))
    MCOLS <- maf[,cols_to_keep,with=FALSE]
    MCOLS[,Genome_Change := paste(Reference_Allele,Tumor_Seq_Allele2,sep='>')]
    MCOLS <- as(MCOLS,'DataFrame')
    maf <- GRanges(seqnames=maf$Chromosome,ranges=IRanges(start=maf$Start_Position,end=maf$End_Position),strand=maf$Strand)
    mcols(maf) <- MCOLS
    if(!is.null(output_file)) {
      saveRDS(maf,file=output_file)
    } else {
      return(maf)
    }
  }
  
  force_seqinfo <- function(x, reference, chromosomes_to_keep=NULL) {
    stopifnot(is(x,'GRanges'))
    si <- seqinfo(reference)
    ss <- seqlevelsStyle(si)
    seqlevelsStyle(x) <- ss
    sl <- seqlevels(si)
    sl <- intersect(sl,seqlevels(x))
    if(!is.null(chromosomes_to_keep)) {
      sl <- intersect(sl,chromosomes_to_keep)
    }
    seqlevels(si) <- sl
    seqlevels(x,pruning.mode='coarse') <- sl
    seqinfo(x) <- si
    return(x)
  }
  
  format_maf_gr <- function(gr, reference, chromosomes_to_keep=NULL) {
    if(inherits(gr,'character')) {
      stopifnot(toupper(tools::file_ext(gr))=='RDS')
      stopifnot(file.exists(gr))
      gr <- readRDS(gr)
    }
    gr <- force_seqinfo(gr,reference=reference,chromosomes_to_keep=chromosomes_to_keep)
    cn <- colnames(mcols(gr))
    cn[grep('project_code',cn,ignore.case=TRUE)] <- 'Project_Code'
    cn[grep('Donor_ID',cn,ignore.case=TRUE)] <- 'Donor_ID'
    cn[grep('Tumor_Sample_Barcode',cn,ignore.case=TRUE)] <- 'Tumor_Sample_Barcode'
    cn[grep('Reference_Allele',cn,ignore.case=TRUE)] <- 'Reference_Allele'
    cn[grep('Tuor_Seq_Allele2',cn,ignore.case=TRUE)] <- 'Tumor_Seq_Allele2'
    cn[grep('Variant_Type',cn,ignore.case=TRUE)] <- 'Variant_Type'
    cn[grep('Project_Code',cn,ignore.case=TRUE)] <- 'Project_Code'
    colnames(mcols(gr)) <- cn
    return(gr)
  }
  
  get_DSR_mutation_stats <- function(DSR_stats, maf_gr=NULL, mappability_track=NULL, normalization_seqs=NULL, blacklist=NULL, reference=NULL, chromosomes_to_keep=NULL, chunk_size=20e3, verbose=TRUE) {
    # DSRs
    if(is.character(DSR_stats)) {
      stopifnot('DSR_stats must be a .csv file (or GRanges object)'=tolower(tools::file_ext(DSR_stats))=='csv')
      dsr <- read.csv(DSR_stats)
      all_dsr <- GRanges(dsr[,1])
      mcols(all_dsr) <- DataFrame(dsr[,-1])
    } else {
      stopifnot('DSR_stats must be a GRanges object (or PATH to .csv)'=is(DSR_stats,'GRanges'))
      all_dsr <- DSR_stats
    }
    if(!is.null(reference)) all_dsr <- force_seqinfo(all_dsr,reference=reference,chromosomes_to_keep=chromosomes_to_keep)
    # Normalization seqs
    if(!is.null(normalization_seqs)) {
      stopifnot('If normalization_seqs is used, reference must be a BSgenome object (and not NULL)'=is(reference,'BSgenome'))
      stopifnot('If normalization_seqs is used, mappability_track must be supplied (not NULL)' = !is.null(mappability_track))
    }
    # Mappability
    if(!is.null(mappability_track)) {
      if(is.character(mappability_track)) {
        if(verbose) message('Importing mappability track: ',basename(mappability_track))
        mt <- import(mappability_track)
      } else {
        stopifnot('mappability_track must be a GRanges object (or PATH to mappability file)'=is(mappability_track,'GRanges'))
        mt <- mappability_track
      }
      if(!is.null(reference)) mt <- force_seqinfo(mt,reference=reference,chromosomes_to_keep=chromosomes_to_keep)
    }
    # Mutations
    if(!is.null(maf_gr)) {
      stopifnot(is(maf_gr,'GRanges'))
      if(!is.null(reference)) maf_gr <- force_seqinfo(maf_gr,reference=reference,chromosomes_to_keep=chromosomes_to_keep)
    }
    # Blacklist
    if(!is.null(blacklist)) {
      stopifnot('If blacklist is supplied, must also supply maf_gr' = !is.null(maf_gr))
      if(is.character(blacklist)) {
        stopifnot(file.exists(blacklist))
        bl <- import(blacklist)
      } else {
        stopifnot(is(blacklist,'GRanges'))
        bl <- blacklist
      }
      if(!is.null(reference)) bl <- force_seqinfo(bl,reference=reference,chromosomes_to_keep=chromosomes_to_keep)
      if(verbose) message('Filtering out mutations overlapping blacklisted regions in ',blacklist)
      idx <- overlapsAny(maf_gr,bl)
      maf_gr <- maf_gr[!idx]
      if(verbose) message('Removed ',round(100*(sum(idx)/length(idx)),digits=1),'% of mutations due to blacklist')
    }
    # Process
    if(verbose) message('Processing ...')
    if(!is.null(maf_gr)) {
      if(verbose) message('Counting mutations in DSR windows ...')
      if(any(colnames(mcols(all_dsr))=='mutation_count')) {
        if(verbose) message('column "mutation_count" already exists and will be overwritten')
        mcols(all_dsr)['mutation_count'] <- NULL
      }
      all_dsr <- binnedSum(all_dsr,coverage(maf_gr,weight=1),'mutation_count')
    }
    if(!is.null(mappability_track)) {
      #if(verbose) message('Counting mappable bases in DSR windows ...')
      all_dsr <- mappable_total(all_dsr,mappability_track=mt,verbose=verbose)
      if(!is.null(normalization_seqs)) {
        if(!is.null(reference)) {
          GEN <- reference
        } else {
          GEN <- all_dsr
        }
        #if(verbose) message('Counting ',paste(normalization_seqs, collapse=', '),' base(s) that overlap mappability')
        all_dsr <- mappable_sequence_counts(all_dsr,
                                            genome=GEN,
                                            mappability_track=mt,
                                            sequences=normalization_seqs,
                                            chunk_size=chunk_size,verbose=verbose)
      }
    }
    if(verbose) message('Done')
    return(all_dsr)
  }
  
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
  
  get_DSR_mutation_stats_by_tumor <- function(DSR_stats, maf_gr, mappability_track=NULL, normalization_seqs=NULL, blacklist=NULL, reference=NULL, chromosomes_to_keep=NULL, chunk_size=20e3, verbose=TRUE, nproc=1) {
    stopifnot('maf_gr must be a GRanges object in maf_gr format'=is(maf_gr,'GRanges'))
    tumors <- as.character(unique(maf_gr$Tumor_Sample_Barcode))
    if(verbose) message('Processing mutations from ',length(tumors),' unique tumors')
    if(!is.null(mappability_track)) {
      DSR_stats <- get_DSR_mutation_stats(DSR_stats=DSR_stats,
                                          maf_gr=NULL,
                                          mappability_track=mappability_track,
                                          normalization_seqs=normalization_seqs,
                                          blacklist=NULL,
                                          reference=reference,
                                          chromosomes_to_keep=chromosomes_to_keep,
                                          chunk_size=chunk_size,
                                          verbose=verbose)
    }
    dat_tumor <- mclapply(tumors,function(i) {
      if(verbose) message('Processing ',i,'     \r',appendLF=FALSE)
      flush.console()
      maf_gri <- maf_gr[maf_gr$Tumor_Sample_Barcode==i]
      get_DSR_mutation_stats(DSR_stats=DSR_stats,
                             maf_gr=maf_gri,
                             mappability_track=NULL,
                             normalization_seqs=NULL,
                             blacklist=blacklist,
                             reference=reference,
                             chromosomes_to_keep=chromosomes_to_keep,
                             chunk_size=chunk_size,
                             verbose=verbose)
    },mc.cores=nproc)
    names(dat_tumor) <- tumors
    dat_tumor <- GRangesList(dat_tumor,compress=TRUE)
    if(verbose) message('\nDONE.')
    return(dat_tumor)
  }
  
  
  
  
  
  ### Define key files
  # DSR statistics
  DSR_stats <- '~/Desktop/trad/DSR_TP53_50kb_resASH.csv' # Generated in this study
  # MAF files
  tcga_maf <-     'MAF/final_consensus_passonly.snv_mnv_indel.tcga.controlled.maf.gz' # downloaded from https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel
  icgc_maf <-     'MAF/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz' # downloaded from https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel
  # SSM files
  skca_br_ssm <-  'SSM/simple_somatic_mutation.open.SKCA-BR.tsv.gz' # downloaded from https://dcc.icgc.org/api/v1/download?fn=/current/Projects/SKCA-BR/simple_somatic_mutation.open.SKCA-BR.tsv.gz
  mela_au_ssm <-  'SSM/simple_somatic_mutation.open.MELA-AU.tsv.gz' # downloaded from https://dcc.icgc.org/api/v1/download?fn=/current/Projects/MELA-AU/simple_somatic_mutation.open.MELA-AU.tsv.gz
  # Miscellany
  mappability <-  'MISC/wgEncodeCrgMapabilityAlign100mer.bigWig' # downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig
  blacklist <-    'MISC/hg19-blacklist.v2.bed.gz' # downloaded from https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz
  # Genome Reference
  REF <- BSgenome.Hsapiens.UCSC.hg19
  chromosomes_to_keep <- seqlevels(REF)[1:23]
  
  
  
  ### Prepare data
  # Convert SSM files to MAF
  maftools::icgcSimpleMutationToMAF(icgc=skca_br_ssm,basename='MAF/SKCA-BR',removeDuplicatedVariants=TRUE)
  maftools::icgcSimpleMutationToMAF(icgc=mela_au_ssm,basename='MAF/MELA-AU',removeDuplicatedVariants=TRUE)
  system('gzip MAF/SKCA-BR.maf MAF/MELA-AU.maf')
  # Convert MAF files to GenomicRanges
  convert_maf_to_granges(maf=tcga_maf,output_file='RDS/TCGA_GRanges.RDS')
  convert_maf_to_granges(maf='MAF/SKCA-BR.maf.gz',output_file='RDS/SKCA-BR_GRanges.RDS')
  convert_maf_to_granges(maf='MAF/MELA-AU.maf.gz',output_file='RDS/MELA-AU_GRanges.RDS')
  # Get Mutational Signatures
  mutation_signatures_from_maf(tcga_maf,output_file='MSIG/TCGA_msig.csv')
  mutation_signatures_from_maf('MAF/SKCA-BR.maf.gz',output_file='MSIG/SKCA-BR_msig.csv')
  mutation_signatures_from_maf('MAF/MELA-AU.maf.gz',output_file='MSIG/MELA-AU_msig.csv')
  
  
  
  ### Import data
  tcga <- format_maf_gr('RDS/TCGA_GRanges.RDS',reference=REF,chromosomes_to_keep=chromosomes_to_keep)
  skcm <- tcga[tcga$Project_Code=='Skin-Melanoma']
  skca <- format_maf_gr('RDS/SKCA-BR_GRanges.RDS',reference=REF,chromosomes_to_keep=chromosomes_to_keep)
  mela <- format_maf_gr('RDS/MELA-AU_GRanges.RDS',reference=REF,chromosomes_to_keep=chromosomes_to_keep)
  tcga.msig <- read.csv('MSIG/TCGA_msig.csv')
  skcm.msig <- tcga.msig[tcga.msig$X %in% unique(skcm$Tumor_Sample_Barcode),]
  skca.msig <- read.csv('MSIG/SKCA-BR_msig.csv')
  mela.msig <- read.csv('MSIG/MELA-AU_msig.csv')
  
  

  ### Define GRanges object for processing 
  ### NOTE, datasets can be combined and analyzed or analyzed individually
  ### Combined
  maf_gr <- c(skcm,skca,mela)
  msig <- rbind(skcm.msig,skca.msig,mela.msig)
  ### OR Individual
  maf_gr <- skcm
  msig <- skcm.msig
  

  
  ### Process
  # Single nucleotide changes only
  maf_gr <- maf_gr[!is.na(maf_gr$Variant_Type) & maf_gr$Variant_Type=='SNP']
  # Sanity check (make sure correct genome assembly)
  seqs <- getSeq(REF,maf_gr)
  stopifnot(all(seqs == maf_gr$Reference_Allele))
  # Sanity check (tumor barcodes identical between file types)
  stopifnot(all(msig$X %in% unique(maf_gr$Tumor_Sample_Barcode)))
  # Filter on C>T mutations
  maf_gr <- maf_gr[maf_gr$Genome_Change %in% c('C>T','G>A')]
  # Drop unused factor levels
  maf_gr$Tumor_Sample_Barcode <- droplevels(maf_gr$Tumor_Sample_Barcode)
  # Import DSR stats
  dsr <- read.csv(DSR_stats)
  MCOLS <- dsr[,-1]
  dsr <- GRanges(dsr[,1])
  mcols(dsr) <- MCOLS
  dsr <- force_seqinfo(dsr,reference=REF,chromosomes_to_keep=chromosomes_to_keep)
  # Calculate DSR mutation stats (total number of mappable bases, number of Cs/Gs overlapping mappability, number of mutations)
  dsr <- get_DSR_mutation_stats(dsr,
                                maf_gr=maf_gr, # C>T mutations in MAF->GRanges format
                                mappability_track=mappability, # Mappability track
                                normalization_seqs=c('C','G'), # Sequences to tally per bin (overlapping mappability)
                                blacklist=blacklist, # Blacklisted regions for excluding mutation counts
                                reference=REF, # BSgenome object (genome reference)
                                chromosomes_to_keep=chromosomes_to_keep) # chromosomes to include for analysis
  # Filter uninformative windows (0 mutations AND 0 CPDs)
  dsr <- dsr[dsr$baseMean > 0 & dsr$mutation_count > 0]
  # Calculate mutation counts per tumor
  tumors <- get_DSR_mutation_stats_by_tumor(dsr,
                                            maf_gr=maf_gr,
                                            mappability_track=NULL, # NULL here since previously computed
                                            normalization_seqs=NULL, # NULL here since previously computed
                                            blacklist=blacklist,
                                            reference=REF,
                                            chromosomes_to_keep=chromosomes_to_keep,
                                            nproc=1)
  # STOPPING POINT
  # dsr <- readRDS('RDS/dsr_03252023.RDS')
  # maf_gr <- readRDS('RDS/maf_gr_03262023.RDS')
  # tumors <- readRDS('RDS/tumors_03252023.RDS')
  # msig <- readRDS('RDS/msig_03252023.RDS')
  # Calculate mutation rate (Number of C>T mutations per 1000 mappable C/Gs) and additional metrics
  mr <- function(dsr) {
    mcols(dsr)$mutation_rate <- (dsr$mutation_count / (dsr$mappable_C + dsr$mappable_G))*500 
    dsr
  }
  dsr <- mr(dsr)
  tumors <- endoapply(tumors,mr)
  # sanity check
  stopifnot(identical(granges(dsr),granges(tumors[[1]])))
  idx <- dsr$padj < 0.0001 & dsr$log2FoldChange > 1
  # solar signature
  uv <- setNames(rowSums(msig[,c('SBS7a','SBS7b','SBS7c','SBS7d')]) / rowSums(msig[,-1]),msig[,1])
  hist(uv,breaks=500)
  uv.bin <- data.frame(Solar_Signature=ifelse(uv > 0.5,'UV','non-UV'))
  uv.bin$Solar_Signature <- factor(uv.bin$Solar_Signature,levels=c('UV','non-UV'))
  tumors.UV <- row.names(uv.bin[uv.bin$Solar_Signature=='UV',,drop=FALSE])
  tumors.nonUV <- row.names(uv.bin[uv.bin$Solar_Signature=='non-UV',,drop=FALSE])
  # percent change
  pc <- function(dsr,idx) (( mean(dsr[idx]$mutation_rate) - mean(dsr$mutation_rate)) / mean(dsr$mutation_rate))*100
  percent_change <- sapply(tumors,pc,idx)
  df2 <- merge(uv.bin,data.frame(percent_change=percent_change),by='row.names')
  ggplot(df2,aes(x=1,y=percent_change)) +
    geom_violin() +
    geom_jitter(aes(color=Solar_Signature),shape=16, position=position_jitter(0.4)) +
    scale_color_brewer(palette='Dark2',direction=-1) +
    # geom_hline(yintercept=0,linetype='dashed') +
    theme_bw()
  # wilcox.test
  wcx <- function(dsr,idx) wilcox.test(dsr[idx]$mutation_rate,dsr[!idx]$mutation_rate)$p.value
  WCX <- sapply(tumors,wcx,idx)
  WCX <- p.adjust(WCX,method='BH')
  stopifnot(identical(names(percent_change),names(WCX)))
  WCX <- sign(percent_change) * -1*log10(WCX)
  df <- merge(uv.bin,data.frame(Pvalue=WCX),by='row.names')
  ggplot(df,aes(x=1,y=Pvalue)) +
    geom_violin() +
    geom_jitter(aes(color=Solar_Signature),shape=16, position=position_jitter(0.4)) +
    scale_color_brewer(palette='Dark2',direction=-1) +
    theme_bw()
  # Mean mutation rates
  gmr <- sapply(tumors,function(i) mean(i[!idx]$mutation_rate)) # rest of genome
  dmr <- sapply(tumors,function(i) mean(i[idx]$mutation_rate))
  # Mutation rate in regions near DSRs
  DSR_width <- median(width(dsr))
  bin_size <- 100
  flanking <- 250e3
  n_bins_dsr <- DSR_width / bin_size
  n_bins_flanking <- flanking / bin_size
  profile <- aggregate_bigwig(dsr[idx],bigwig_file=coverage(maf_gr),bin_size=bin_size,flanking=flanking,ignore.strand=TRUE)
  plot(profile,pch=19,col=rgb(t(col2rgb('gray50')),maxColorValue=255,alpha=100))
  lines(smooth.spline(profile,spar=0.7),col='darkred',lwd=3,lty='dashed')
  abline(v=c(0,n_bins_flanking,n_bins_flanking+n_bins_dsr,n_bins_flanking*2+n_bins_dsr),lty='dotted')
  # Dot plot
  #order <- names(dmr[order(dmr)])
  #order <- names(gmr[order(gmr)])
  
  WCX <- WCX[order(WCX)]
  order <- names(WCX)
  ptl_line <- which(WCX > -log10(0.05))[1] # partition line
  stopifnot(identical(names(gmr),names(dmr)))
  df3 <- data.frame(Tumor=names(gmr),Rest_of_Genome=gmr,DSR=dmr)
  df3$Tumor <- factor(df3$Tumor,levels=order)
  DOT <- pivot_longer(df3,cols=c(Rest_of_Genome,DSR),values_to='Mutation_Rate',names_to='Category')
  ggplot(DOT,aes(x=Tumor,y=Mutation_Rate)) +
    geom_line(aes(group=Tumor),linewidth=0.25) +
    geom_point(aes(color=Category),size=1) +
    geom_vline(xintercept = ptl_line,linetype='dotted') +
    scale_y_continuous(trans='log') +
    scale_color_manual(values=c('Rest_of_Genome'='#A1D76A','DSR'='#E9A3C9')) + 
    scale_x_discrete(labels=NULL,breaks=NULL) + labs(x=NULL) +
    theme_bw()
  pheatmap(WCX[order],cluster_rows=F,cluster_cols=F)
  
  # Mutation signatures
  # top signatures
  ts.uv <- colSums(msig[msig$X %in% tumors.UV,-1])
  ts.nonUV <- colSums(msig[msig$X %in% tumors.nonUV,-1])
  ts.uv <- sort(ts.uv,TRUE)
  ts.nonUV <- sort(ts.nonUV,TRUE)
  barplot(ts.uv,las=2)
  barplot(ts.nonUV,las=2)
  top_sigs <- union(names(ts.uv)[1:8],names(ts.nonUV)[1:8])
  ms <- msig
  row.names(ms) <- ms$X
  ms <- as.data.frame(t(ms[,-1]))
  ms$SBS <- row.names(ms)
  ms <- pivot_longer(ms,cols=!SBS,names_to='Tumor')
  ms$Tumor <- factor(ms$Tumor,levels=order)
  
  uv_sig <- setNames('#FAF39D','#FFCC66','#F58022','red',c('SBS7a','SBS7b','SBS7c','SBS7d'))
  unknown <- setNames(brewer.pal(3,'Greens'),c('SBS43','SBS39','SBS37'))
  seq_art <- setNames(brewer.pal(2,'Purples')[1:2],c('SBS58','SBS50'))
  apobec <- setNames(brewer.pal(1,'Blues')[1],c('SBS2'))
  alkyl <- setNames('pink',c('SBS11'))
  clock <- setNames(brewer.pal(1,'Set3')[1],c('SBS5'))
  endog <- setNames(brewer.pal(1,'Set1')[1],c('SBS1'))
  bxr <- setNames('dodgerblue',c('SBS30'))
  
  Cols <- c(uv_sig,unknown,seq_art,apobec,alkyl,clock,endog,bxr)
  
  
  
  
  ggplot(ms,aes(fill=SBS,y=value,x=Tumor)) +
    geom_bar(position='fill',stat='identity',width=1) +
    #scale_fill_d3(alpha=0.4) +
    scale_fill_manual(values=Cols)
  

    scale_fill_manual(values=c('SBS1'='#A6CEE3','SBS2'='#6A3E98','SBS5'='#C9B3D6','SBS6'='#1A79B4','SBS7a'='#FAF39D','SBS7b'='#FFCC66','SBS7c'='#F58022','SBS7d'='red','SBS8'='#B69064','SBS13'='#B2D88B','SBS18'='darkred'))
  
  
  
  
  
  
  
  
  
  
  # Mutation rate in regions near DSRs UV tumors
  DSR_width <- median(width(dsr))
  bin_size <- 100
  flanking <- 250e3
  n_bins_dsr <- DSR_width / bin_size
  n_bins_flanking <- flanking / bin_size
  profile <- aggregate_bigwig(dsr[idx],bigwig_file=coverage(maf_gr[maf_gr$Tumor_Sample_Barcode %in% tumors.UV]),bin_size=bin_size,flanking=flanking,ignore.strand=TRUE)
  plot(profile,pch=19,col=rgb(t(col2rgb('gray50')),maxColorValue=255,alpha=100))
  lines(smooth.spline(profile,spar=0.7),col='darkred',lwd=3,lty='dashed')
  abline(v=c(0,n_bins_flanking,n_bins_flanking+n_bins_dsr,n_bins_flanking*2+n_bins_dsr),lty='dotted')
  
  # Mutation rate in regions near DSRs UV tumors
  DSR_width <- median(width(dsr))
  bin_size <- 100
  flanking <- 250e3
  n_bins_dsr <- DSR_width / bin_size
  n_bins_flanking <- flanking / bin_size
  profile <- aggregate_bigwig(dsr[idx],bigwig_file=coverage(maf_gr[maf_gr$Tumor_Sample_Barcode %in% tumors.nonUV]),bin_size=bin_size,flanking=flanking,ignore.strand=TRUE)
  plot(profile,pch=19,col=rgb(t(col2rgb('gray50')),maxColorValue=255,alpha=100))
  lines(smooth.spline(profile,spar=0.7),col='darkred',lwd=3,lty='dashed')
  abline(v=c(0,n_bins_flanking,n_bins_flanking+n_bins_dsr,n_bins_flanking*2+n_bins_dsr),lty='dotted')
  
  
  
  #pheatmap(pc[order],cluster_rows=F,cluster_cols=F)
  

  
  
  ggplot(DOT,aes(x=id,y=value)) +
    geom_line(aes(group=id))
  
  ggplot(DOT,aes(value,id)) +
    geom_line(aes(group=id)) +
    geom_point(aes(color=variable),size=2) +
    coord_flip() +
    scale_color_manual(values=c('Genome'='#A1D76A','DSR'='#E9A3C9')) + 
    theme_bw()
  pheatmap(t(PERCENT_CHANGE),cluster_rows=FALSE,cluster_cols=FALSE,color=colorRampPalette(rev(brewer.pal(7,'PiYG')))(100))
  

  
  
  
  # Detailed plots
  pss <- fread(pcawg_sample_sheet)
  ms <- fread(mutational_signatures)
  cn <- fread(copy_number)
  cn <- as.data.frame(cn)
  row.names(cn) <- cn$`Gene Symbol`
  # Add sample barcode to mutational signatures data
  idx <- setNames(pss$submitter_sample_id,pss$icgc_specimen_id)
  ms$Tumor_Sample_Barcode <- idx[ms$`Sample Names`]
  # sanity check
  sum(tumors %in% pss$submitter_sample_id)
  sum(tumors %in% colnames(cn))
  sum(tumors %in% ms$Tumor_Sample_Barcode)
  # p53 pathway genes
  p53_path_genes <- c('CHEK2','CHEK1','ATR','ATM','TP53','MDM4','MDM2')
  all(p53_path_genes %in% cn$`Gene Symbol`)
  p53_path_gene_IDs <- select(org.Hs.eg.db,columns='ENTREZID',keys=p53_path_genes,keytype='SYMBOL')
  gene_locs <- genes(Homo.sapiens)
  gene_locs <- gene_locs[unlist(gene_locs$GENEID %in% as.character(p53_path_gene_IDs$ENTREZID))]
  gene_locs$SYMBOL <- mapIds(Homo.sapiens,keys=unlist(gene_locs$GENEID),column='SYMBOL',keytype='ENTREZID')
  gene_locs <- keepSeqlevels(gene_locs,value=seqlevels(gene_locs)[1:23],pruning.mode='coarse')
  ms_codes <- c('SBS1','SBS2','SBS5','SBS6','SBS7a','SBS7b','SBS7c','SBS7d','SBS8','SBS13','SBS18')
  
  dat <- lapply(tumors,function(i) {
    message(i)
    tumor_all <- mut[mut$Tumor_Sample_Barcode==i]
    tumor_ctot <- ctot[ctot$Tumor_Sample_Barcode==i]
    tumor_ms <- as.data.frame(ms[ms$Tumor_Sample_Barcode==i])
    tumor_cn <- t(cn[p53_path_genes,colnames(cn)==i,drop=FALSE])
    if(length(tumor_ctot) < 100) {
      list(GEN_CTOT=NA,DSR_CTOT=NA,PERCENT_CHANGE=NA,CTOT_FRAC=NA,P53_PATH_SSM=NA,TUMOR_MS=NA)
    } else {
      # C>T mutation rate in genome bins and top DSRs
          geni <- all_dsr
          geni <- binnedAverage(geni,coverage(tumor_ctot),'mutations.i')
          geni$mutations.i <- geni$mutations.i*width(geni)
          geni$mutation_rate.i <- geni$mutations.i / (geni$mappable_C + geni$mappable_G)
          dsri <- geni[geni$padj < padj_threshold & geni$log2FoldChange > lfc_threshold]
          GEN_CTOT <- mean(geni$mutation_rate.i)*1e4
          DSR_CTOT <- mean(dsri$mutation_rate.i)*1e4
          PERCENT_CHANGE <- ((DSR_CTOT - GEN_CTOT)/GEN_CTOT)*100
      # C>T fraction
          CTOT_FRAC <- length(tumor_ctot) / length(tumor_all)
      # p53 pathway SSM
          ssm_genes <- binnedSum(gene_locs,numvar=coverage(tumor_all),varname='x')
          P53_PATH_SSM <- setNames(ssm_genes$x,nm=ssm_genes$SYMBOL)
      # Mutatational Signature
          TUMOR_MS <- tumor_ms[,colnames(tumor_ms) %in% ms_codes]
      # Copy Number
          TUMOR_CN <- as.data.frame(tumor_cn)
      # Return results
          list(GEN_CTOT=GEN_CTOT,DSR_CTOT=DSR_CTOT,PERCENT_CHANGE=PERCENT_CHANGE,CTOT_FRAC=CTOT_FRAC,P53_PATH_SSM=P53_PATH_SSM,TUMOR_MS=TUMOR_MS,TUMOR_CN=TUMOR_CN)
    }
  })
  names(dat) <- tumors
  ORDER <- order(sapply(dat,'[[','DSR_CTOT'))
  GEN_CTOT <- sapply(dat,'[[','GEN_CTOT')[ORDER]
  DSR_CTOT <- sapply(dat,'[[','DSR_CTOT')[ORDER]
  PERCENT_CHANGE <- sapply(dat,'[[','PERCENT_CHANGE')[ORDER]
  CTOT_FRAC <- sapply(dat,'[[','CTOT_FRAC')[ORDER]
  P53_PATH_SSM <- sapply(dat,'[[','P53_PATH_SSM')[,ORDER]
  TUMOR_MS <- sapply(dat,'[[','TUMOR_MS')[,ORDER]
  TUMOR_CN <- sapply(dat,'[[','TUMOR_CN')[,ORDER]
  TUMOR_CN <- apply(TUMOR_CN,2,unlist)
  MS <- as.data.frame(t(TUMOR_MS))
  MS$id <- row.names(MS)
  MS$id <- factor(MS$id,levels=tumors[ORDER])
  
  MS <- as.data.frame(MS %>% 
    pivot_longer(ms_codes))
  MS$value <- unlist(MS$value)
  
  DOT <- data.frame(value=c(unname(GEN_CTOT),unname(DSR_CTOT)),id=c(names(GEN_CTOT),names(DSR_CTOT)),variable=c(rep('Genome',37),rep('TP53 DSR',37)))
  DOT$id <- factor(DOT$id,levels=names(GEN_CTOT))
  
  
  ggplot(DOT,aes(value,id)) +
    geom_line(aes(group=id)) +
    geom_point(aes(color=variable),size=2) +
    coord_flip() +
    scale_color_manual(values=c('Genome'='#A1D76A','TP53 DSR'='#E9A3C9')) + 
    theme_bw()
  pheatmap(t(PERCENT_CHANGE),cluster_rows=FALSE,cluster_cols=FALSE,color=colorRampPalette(rev(brewer.pal(7,'PiYG')))(100))
  pheatmap(P53_PATH_SSM,cluster_rows=FALSE,cluster_cols=FALSE,color=colorRampPalette(rev(brewer.pal(n=7,name='RdBu')))(7),
           border_color=NA,
           breaks=seq(-5,5,length.out=8))
  pheatmap(TUMOR_CN,cluster_rows=FALSE,cluster_cols=FALSE,color=colorRampPalette(rev(brewer.pal(n=7,name='BrBG')))(7),
           border_color=NA)
  pheatmap(t(CTOT_FRAC),cluster_rows=FALSE,cluster_cols=FALSE,color=colorRampPalette(brewer.pal(7,'Purples'))(100) )
  ggplot(MS,aes(fill=name,y=value,x=id)) +
    geom_bar(position='fill',stat='identity',width=1) +
    #scale_fill_d3(alpha=0.4)
    scale_fill_manual(values=c('SBS1'='#A6CEE3','SBS2'='#6A3E98','SBS5'='#C9B3D6','SBS6'='#1A79B4','SBS7a'='#FAF39D','SBS7b'='#FFCC66','SBS7c'='#F58022','SBS7d'='red','SBS8'='#B69064','SBS13'='#B2D88B','SBS18'='darkred'))

  
  
  
  
  
  
  
  
  
  ## Dot plot
  df <- data.frame(donor=dat$icgc_donor_id, gen=dat$genome_mean_mutation_rate_50kb_10kb,dsr=dat$TP53DSR_mean_mutation_rate_50kb_10kb)
  #df <- df[order(df$dsr),]
  dm <- reshape2::melt(df,id='donor')
  dm$donor <- factor(dm$donor,levels=donor_order)
  dm$variable <- relevel(dm$variable,ref='dsr')
  ggplot(dm,aes(value,donor)) +
    geom_line(aes(group=donor)) +
    geom_point(aes(color=variable),size=2) +
    coord_flip() +
    scale_color_manual(values=c('gen'='#A1D76A','dsr'='#E9A3C9')) + 
    theme_bw()

  
  
  

  
  
  

 