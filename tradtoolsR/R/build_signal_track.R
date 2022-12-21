
  # input_bam <- '~/Desktop/trad/old/misc/IMR_100J_UVC_1_duprm_spec1.bam'
  # input_blacklist <- '~/Desktop/trad/old/misc/hg19-blacklist.v2.bed.gz'
  # output_bed <- '~/Desktop/trad/old/misc/tmppp.bed'
  # output_bw <- '~/Desktop/trad/old/misc/tmppp.bw'
  # export_stats <- T
  # reference_fasta <- '~/Desktop/trad/old/misc/spec1.fa'
  # ranges_operation <- 'default'
  # mapQ <- 10
  # standard_chromosomes <- F
  # overwrite <- T
  # verbose <- T

  # setwd('~/Desktop/trad/tradtools/tradtools/')
  # devtools::document()
  # setwd('~/Desktop/trad/misc')
  # build_signal_track(input_bam = 'IMR_100J_UVC_1_duprm_spec1.bam', input_blacklist = 'hg19-blacklist.v2.bed.gz', output_bed='tmp.bed', output_bw='tmp.bw',reference_fasta = 'spec1.fa', overwrite=T)


  #' @title TRADtools build_signal_track
  #' @description This function processes a BAM file into a TRAD-seq signal track file in BED format
  #' @usage \code{build_signal_track(...)}
  #' @param input_bam PATH to the input BAM file. BAM file should be pre-processed and sorted
  #' @param input_blacklist PATH to the blacklist BED file to be used for filtering (e.g. "/path/to/hg19-blacklist.v2.bed.gz"), can be NULL (or "NULL" or "null")
  #' @param output_bed PATH to the output BED file (will not be overwritten if already exists, unless overwrite=TRUE)
  #' @param output_bw PATH to the output BW file (will not be overwritten if already exists, unless overwrite =T)
  #' @param export_stats TRUE or FALSE. Whether to write out stats files from the signal track processing (to the current directory)
  #' @param reference_fasta PATH to the reference FASTA file that was used for alignment. This provides the sequences that will be retrieved after applying ranges_operation() to the BAM reads
  #' @param ranges_operation can be "default" (two bases immediately 5' to R1 read), "identity" (no transformation), or a user-defined function
  #' @param mapQ Integer. Reads with mapQ below this number will be discarded
  #' @param standard_chromosomes TRUE or FALSE. If TRUE, retain only the standard chromosomes (by using \code{GenomeInfoDb::keepStandardChromosomes()})
  #' @param overwrite TRUE or FALSE. If TRUE, output BED file will be overwritten even if it already exists
  #' @param verbose TRUE or FALSE. Whether to report helpful messages and summary statistics
  #' @details The exported BED file contains reads after applying mapQ filter and transformation defined by ranges_operation. The exported BW file reflects the "final" signal after further filtering out non-dipyrimidine records and records that overlap a black-listed region.
  #' @return Nothing returned. Files are exported.
  #' @examples
  #' \code{build_signal_track(gen1.bam, bl.bed.gz, gen1.bed, gen1.bw, gen1.fasta)}
  #' @import GenomicRanges
  #' @import BSgenome
  #' @import Biostrings
  #' @import Rsamtools
  #' @import GenomicAlignments
  #' @import GenomeInfoDb
  #' @import rtracklayer
  #' @export
  build_signal_track <- function(input_bam, input_blacklist, output_bed, output_bw, export_stats=TRUE, reference_fasta, ranges_operation='default', mapQ=10, standard_chromosomes=FALSE, overwrite=FALSE, verbose=TRUE) {
    t0 <- Sys.time()
    # ------ check args  ----- #
    if(!file.exists(input_bam)) stop('File "',input_bam,'" not found')
    if(is.null(input_blacklist) | toupper(basename(input_blacklist))=='NULL') {
      input_blacklist <- 'NULL'
    } else {
      if(!file.exists(input_blacklist)) stop('File "',input_blacklist,'" not found')
    }
    if(!file.exists(reference_fasta)) stop('File "',reference_fasta,'" not found')
    if(!is.logical(standard_chromosomes)) stop('standard_chromosomes must be a logical value (TRUE or FALSE)')
    if(!is.logical(overwrite)) stop('overwrite must be a logical value (TRUE or FALSE)')
    if(!is.logical(verbose)) stop('verbose must be a logical value (TRUE or FALSE)')
    if(!is.logical(export_stats)) stop('export_stats must be a logical value (TRUE or FALSE)')
    if(!(is.numeric(mapQ) & length(mapQ)==1)) stop('mapQ must be a single number')
    #if(!(mapQ >= 0 & mapQ <= 42)) stop('mapQ must be between 0 and 42')
    if(overwrite) {
      if(file.exists(output_bed)) {
        if(verbose) message('Overwriting output file ',output_bed)
        unlink(output_bed)
      }
      if(file.exists(output_bw)) {
        if(verbose) message('Overwriting output file ',output_bw)
        unlink(output_bw)
      }
    } else {
      if(file.exists(output_bed)) stop('File "',output_bed,'" already exists and will not be overwritten. Use overwrite = TRUE to force overwriting the output file')
      if(file.exists(output_bw)) stop('File "',output_bw,'" already exists and will not be overwritten. Use overwrite = TRUE to force overwriting the output file')
    }
    if(is.character(ranges_operation)) {
      if(!(ranges_operation %in% c('default','identity'))) stop('If ranges_operation is provided as a character string, it must be one of "default" or "identity"')
      ranges_operation <- ifelse(ranges_operation=='default',yes = function(gr) GenomicRanges::flank(gr,width=2,start=T,both=F,use.names=F,ignore.strand=F), no = function(x)x)
    }
    if(!is.function(ranges_operation)) stop('If ranges_operation is not "default" or "identity", it must be a function')
    test_gr <- GenomicRanges::GRanges('chr1:10000-20000:+')
    tryCatch(
      a <- ranges_operation(test_gr), error=function(e) stop('The supplied ranges_operation function is not valid. See examples for proper usage')
    )
    # ----- import FASTA ----- #
    if(verbose) message('Importing reference fasta file "',basename(reference_fasta),'" ...',appendLF=F)
    ref <- Biostrings::readDNAStringSet(reference_fasta)
    if(verbose) message('\rImporting reference fasta file "',basename(reference_fasta),'" ... DONE.')
    # ----- import BAM   ----- #
    if(verbose) message('Importing BAM file "',basename(input_bam),'" ...',appendLF=F)
    LOG_R1_number_initial <- Rsamtools::countBam(input_bam,param=Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isFirstMateRead=TRUE)))$records
    #LOG_R2_number_initial <- Rsamtools::countBam(input_bam,param=Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isSecondMateRead=TRUE)))$records
    sbf <- Rsamtools::scanBamFlag(isFirstMateRead=T, isProperPair=T, isUnmappedQuery=F, isSecondaryAlignment=F, isNotPassingQualityControls=F, isDuplicate=F) # only need to read R1
    sbp <- Rsamtools::ScanBamParam(mapqFilter=mapQ,flag=sbf)
    loc <- GenomicAlignments::readGAlignments(input_bam, param=sbp)
    if(verbose) message('\rImporting BAM file "',basename(input_bam),'" ... DONE.')
    if(verbose) message('Processing BAM file "',basename(input_bam),'" ...',appendLF=F)
    loc <- as(loc,'GenomicRanges')
    LOG_RN_after_input <- length(loc)
    loc <- ranges_operation(loc)
    invisible(gc())
    #loc <- GenomicRanges::flank(loc,width=2,start=T,both=F,use.names=F,ignore.strand=F) # transform the ranges
    loc <- GenomicRanges::trim(loc)
    wl <- GenomicRanges::width(loc)
    .mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x,ux)))]
    }
    loc <- loc[wl==.mode(wl)]
    rm(wl)
    invisible(gc())
    # --- remove blacklist --- #
    if(input_blacklist!='NULL') {
      bl <- rtracklayer::import.bed(input_blacklist)
      shared_seqinfo2 <- GenomeInfoDb::intersect(GenomeInfoDb::seqinfo(bl),GenomeInfoDb::seqinfo(loc))
      if(length(shared_seqinfo2)==0) stop('There are 0 matching chromosome names in the BAM file and blacklist file. Make sure the chromosome naming conventions are consistent')
      GenomeInfoDb::seqlevels(loc,pruning.mode='coarse') <- GenomeInfoDb::seqlevels(shared_seqinfo2)
      loc <- loc[!IRanges::overlapsAny(loc, bl,ignore.strand=T)] # filter out black-listed records
    }
    # ---- add sequences ----- #
    shared_seqinfo <- GenomeInfoDb::intersect(GenomeInfoDb::seqinfo(ref),GenomeInfoDb::seqinfo(loc))
    if(length(shared_seqinfo)==0) stop('There are 0 matching chromosome names in the BAM file and reference file. Make sure the reference is the same one used for alignment')
    GenomeInfoDb::seqlevels(loc,pruning.mode='coarse') <- GenomeInfoDb::seqlevels(shared_seqinfo)
    if(standard_chromosomes) loc <- GenomeInfoDb::keepStandardChromosomes(loc,pruning.mode='coarse')
    if(length(loc)==0) stop('There are 0 reads remaining')
    seqs <- BSgenome::getSeq(ref,loc)
    seqs <- as.character(Biostrings::reverseComplement(seqs))
    GenomicRanges::mcols(loc)['name'] <- seqs
    GenomicRanges::mcols(loc)['score'] <- 1L
    LOG_RN_after_processing <- length(loc)
    LOG_seqs_after_processing <- table(seqs,dnn=NULL)
    seq_order <- c('AA', 'AG', 'AC', 'AT', 'GA', 'GG', 'GC', 'GT', 'CA', 'CG', 'CC', 'CT', 'TA', 'TG', 'TC', 'TT')
    LOG_seqs_after_processing <- LOG_seqs_after_processing[seq_order]
    names(LOG_seqs_after_processing) <- seq_order
    LOG_seqs_after_processing[is.na(LOG_seqs_after_processing)] <- 0L
    invisible(gc())
    if(verbose) message('\rProcessing BAM file "',basename(input_bam),'" ... DONE.')
    # ------ export BED  ----- #
    if(verbose) message('Exporting BED file "',basename(output_bed),'" ...',appendLF=F)
    rtracklayer::export.bed(object=loc,con=output_bed)
    if(verbose) message('\rExporting BED file "',basename(output_bed),'" ... DONE.')
    # ------ export BW ------- #
    if(verbose) message('Exporting BW file "',basename(output_bw),'" ...',appendLF=F)
    loc <- loc[loc$name %in% c('TT','TC','CT','CC')]
    LOG_RN_after_seqfiltering <- length(loc)
    loc <- GenomicRanges::coverage(loc,weight='score')
    invisible(gc())
    rtracklayer::export.bw(object=loc,con=output_bw)
    if(verbose) message('\rExporting BW file "',basename(output_bw),'" ... DONE.')
    # ----- export logs ------ #
      ID <- tools:::file_path_sans_ext(basename(input_bam))
      # Dinucs
      Line1dn <- c(ID,names(LOG_seqs_after_processing))
      Line1dn <- paste(Line1dn,collapse=',')
      Line2dn <- c(ID,as.numeric(LOG_seqs_after_processing))
      Line2dn <- paste(Line2dn,collapse=',')
      # Read counts
      reads_removed_during_import       <- LOG_R1_number_initial - LOG_RN_after_input
      reads_removed_during_processing   <- LOG_RN_after_input - LOG_RN_after_processing
      reads_removed_during_seqfiltering <- LOG_RN_after_processing - LOG_RN_after_seqfiltering
      final_counts                      <- LOG_RN_after_seqfiltering
      Line1rc <- c(ID,'final_counts','reads_removed_during_dinucleotide_filtering','reads_removed_during_processing','reads_removed_during_import')
      Line1rc <- paste(Line1rc,collapse=',')
      Line2rc <- c(ID,final_counts,reads_removed_during_seqfiltering,reads_removed_during_processing,reads_removed_during_import)
      Line2rc <- paste(Line2rc,collapse=',')
      if(export_stats) {
        writeLines(c(Line1dn,Line2dn),con=paste0('./',ID,'_trad_dinucs.csv'))
        writeLines(c(Line1rc,Line2rc),con=paste0('./',ID,'_trad_readcounts.csv'))
      }
    # ------------------------ #
    if(verbose) {
      message('| ------------------------------- FINAL STATS ------------------------------- |')
      message('|   Input BAM: ',sprintf('%-61s',basename(input_bam)),'  |')
      message('|   Input REF: ',sprintf('%-61s',basename(reference_fasta)),'  |')
      message('|  Output BED: ',sprintf('%-61s',basename(output_bed)),'  |')
      message('|     Initial: ',sprintf('%-19s',format(LOG_R1_number_initial,scientific=F)),' R1 reads                                   |')
      filter_reads_perc <- sprintf('%3s',round((LOG_RN_after_processing/LOG_R1_number_initial)*100,digits=1))
      final_reads_perc <- sprintf('%3s',round((LOG_RN_after_seqfiltering/LOG_R1_number_initial)*100,digits=1))
      message('|Intermediate: ',sprintf('%-19s',format(LOG_RN_after_processing,scientific=F)),' R1 reads after filtering (',filter_reads_perc,'%)           |')
      message('|       Final: ',sprintf('%-19s',format(LOG_RN_after_seqfiltering,scientific=F)),' R1 reads after filtering (',final_reads_perc,'%)           |')
      TT_reads <- LOG_seqs_after_processing[['TT']]
      TC_reads <- LOG_seqs_after_processing[['TC']]
      CT_reads <- LOG_seqs_after_processing[['CT']]
      CC_reads <- LOG_seqs_after_processing[['CC']]
      message('|          TT: ',sprintf('%-19s',format(TT_reads,scientific=F)),sprintf('%45s','|'))
      message('|          TC: ',sprintf('%-19s',format(TC_reads,scientific=F)),sprintf('%45s','|'))
      message('|          CT: ',sprintf('%-19s',format(CT_reads,scientific=F)),sprintf('%45s','|'))
      message('|          CC: ',sprintf('%-19s',format(CC_reads,scientific=F)),sprintf('%45s','|'))
      signal_ratio <- sprintf('%3s',round((sum(TT_reads,TC_reads,CT_reads,CC_reads)/LOG_RN_after_processing)*100,digits=1))
      message('|      Signal: ',signal_ratio,'%',sprintf('%59s','|'))
      message('| ',sprintf('%75s',paste('Total elapse time: ',format(round(Sys.time()-t0,digits=2)))),' |')
      message('| --------------------------------------------------------------------------- |')
    }
  }


