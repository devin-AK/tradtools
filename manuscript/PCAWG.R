
  setwd('~/Desktop/trad/PCAWG/manuscript/')
  
  
  library(R.utils)
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
  library(Homo.sapiens)
  library(pheatmap)
  library(RColorBrewer)
  

  
  
  # Define important files and parameters
  DSR_stats <- '~/Desktop/trad/DSR_TP53_50kb_resASH.csv' # Generated in this study
  mappability_track <- 'wgEncodeCrgMapabilityAlign100mer.bigWig' # downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig
  pcawg_sample_sheet <- 'pcawg_sample_sheet.tsv' # downloaded from https://dcc.icgc.org/releases/PCAWG/donors_and_biospecimens
  mutational_signatures <- 'PCAWG_sigProfiler_SBS_signatures_in_samples.csv' # downloaded from https://dcc.icgc.org/releases/PCAWG/mutational_signatures/Signatures_in_Samples/SP_Signatures_in_Samples
  tcga_maf <- 'final_consensus_passonly.snv_mnv_indel.tcga.controlled.maf.gz' # downloaded from https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel
  icgc_maf <- 'final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz' # downloaded from https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel
  copy_number <- 'all_samples.consensus_level_calls.by_gene.170214.txt' # downloaded from https://dcc.icgc.org/releases/PCAWG/consensus_cnv/gene_level_calls
  
  
  
  
  # Top DSR parameters
  padj_threshold <- 0.00001
  lfc_threshold <- 1
  
  
  
  
  
  # Import and process melanoma mutations from TCGA
  maf <- fread(tcga_maf)
  mut <- GRanges(seqnames=maf$Chromosome,ranges=IRanges(start=maf$Start_position,end=maf$End_position),strand=maf$Strand,
                 Project_Code=maf$Project_Code,
                 Donor_ID=maf$Donor_ID,
                 Tumor_Sample_Barcode=maf$Tumor_Sample_Barcode,
                 Genome_Change=maf$Genome_Change,
                 Reference_Allele=maf$Reference_Allele,
                 Tumor_Seq_Allele2=maf$Tumor_Seq_Allele2,
                 Variant_Type=maf$Variant_Type)
  seqlevelsStyle(mut) <- 'UCSC'
  mut <- sortSeqlevels(mut)
  seqinfo(mut) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
  mut <- keepSeqlevels(mut,seqlevels(mut)[1:23],pruning.mode='coarse')
  #saveRDS(mut,file='TCGA_maf.RDS')
  rm(maf); gc()
  #mut <- readRDS('TCGA_maf.RDS')
  
  mela <- mut[mut$Project_Code=='Skin-Melanoma' & mut$Variant_Type=='SNP']
  mela$Mutation <- substr(mela$Genome_Change,nchar(mela$Genome_Change) - 2,nchar(mela$Genome_Change))
  # sanity check (make sure correct genome assembly)
  seqs <- getSeq(Hsapiens,mela)
  all(as.character(seqs) == substr(mela$mutation,start = 1,stop=1))
  # Filter on C>T mutations
  ctot <- mela[mela$Mutation %in% c('C>T','G>A')]
  
  
  
  
  # Import and process mappability track
  # download from: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig
  mt <- import(mappability_track)
  mt <- sortSeqlevels(mt)
  seqinfo(mt) <- seqinfo(Hsapiens)
  mt <- keepSeqlevels(mt,seqlevels(mt)[1:23],pruning.mode='coarse')
  
  
  
  # Import and process DSR file
  dsr <- read.csv(DSR_stats)
  all_dsr <- GRanges(dsr$X,seqinfo=seqinfo(mela),padj=dsr$padj,log2FoldChange=dsr$log2FoldChange)
  all_dsr <- mappable_sequence_counts(gr=all_dsr,
                                      genome=Hsapiens,
                                      mappability_track=mt,
                                      sequences=c('C','G'),
                                      chunk_size=20e3)
  all_dsr <- mappable_total(all_dsr,mappability_track=mt)
  all_dsr <- binnedAverage(all_dsr,coverage(ctot,weight=1),'mutations')
  all_dsr$mutations <- all_dsr$mutations * width(all_dsr)
  all_dsr$mutation_rate <- all_dsr$mutations / (all_dsr$mappable_C + all_dsr$mappable_G)
  all_dsr <- all_dsr[all_dsr$mutations>0]
  #saveRDS(all_dsr,file='all_dsr.RDS')
  #all_dsr <- readRDS('all_dsr.RDS')
  top_dsr <- all_dsr[all_dsr$padj < padj_threshold & all_dsr$log2FoldChange > lfc_threshold]
  
  
  # C>T mutation rate in DSRs vs genome
  L <- list(Genome=all_dsr$mutation_rate*1e4,DSRs=top_dsr$mutation_rate*1e4)
  cool_vioplot(L) + scale_y_log10() + theme(legend.position = 'none')
  
  
  # C>T mutation rate around DSRs
  cp <- composite_plot(top_dsr, var=coverage(ctot), num_bins=1000, upstream=500e3, downstream=500e3, col1='gray80', col2='darkred', spar=0.7) 
 
  
  # Percent change of TP53-DSR C>T mutation rate vs tumor average in melanoma
  tumors <- unique(ctot$Tumor_Sample_Barcode)
  pc <- lapply(tumors,function(i){
    ctoti <- ctot[ctot$Tumor_Sample_Barcode==i]
    if(length(ctoti)<100) {
      NA
    } else {
      geni <- all_dsr
      geni <- binnedAverage(geni,coverage(ctoti),'mutations.i')
      geni$mutations.i <- geni$mutations.i*width(geni)
      geni$mutation_rate.i <- geni$mutations.i / (geni$mappable_C + geni$mappable_G)
      dsri <- geni[geni$padj < padj_threshold & geni$log2FoldChange > lfc_threshold]
      gmri <- mean(geni$mutation_rate.i,na.rm=TRUE)
      dmri <- mean(dsri$mutation_rate.i,na.rm=TRUE)
      ((dmri - gmri)/gmri)*100
    }
  })
  names(pc) <- tumors
  cool_vioplot(unlist(pc)) + 
    geom_point(position=position_jitter(width=0.3,height=0),size=0.8) + 
    geom_hline(yintercept=0,linetype='dashed') + 
    theme(legend.position='none') + 
    scale_size(range=c(0.1,0.1)) +
    scale_fill_manual(values=add_transparency('#E5F5E0',0.1))
  
  
  
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

  
  
  

  
  
  
  
  
  ### Functions
  
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
    colnames(nuc_counts) <- paste0(prefix,colnames(nuc_counts))
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
    if(verbose) message('Calculating total number of mappable bases in granges ...',appendLF=F)
    mtc <- coverage(mappability_track)
    mtc <- mtc[seqlevels(gr)]
    gr$mappable_total <- GenomicRanges::binnedAverage(bins=gr,numvar=mtc,'x')$x * width(gr)
    message(' DONE.')
    return(gr)
  }
  
  
  cool_vioplot <- function(..., include_points=F, vio_alpha=1, alpha=0.3, main=NULL, xlab=NULL, ylab=NULL, palette='Blues') {
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
      p <- p + geom_jitter(shape=16, position=position_jitter(0.4), alpha=alpha)
    }
    p + geom_boxplot(width=0.1,fill='white',outlier.shape=NA)
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
  
  composite_plot <- function(x, var, num_bins=NULL, bin_size=NULL, upstream=0, downstream=0, ignore.strand=F, spar=NULL, col1='skyblue1', col2='tan1', verbose=TRUE, xlab='', ylab='', ...) {
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
    plot(res,type='l',xaxt='n',ylab=ylab,xlab=xlab,col=col1,lwd=0.25,...)
    lines(smooth.spline(res,spar=spar),col=col2,lwd=2)
    axis(1,at=at,labels=lab,las=2)
    abline(v=c(length(upstr),length(upstr)+length(center)),lty=2)
    
    return(res)
    
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
  
  
 