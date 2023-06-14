
  setwd('~/Desktop/trad')
  
  library(rtracklayer)
  library(BRGenomics)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  library(dplyr)
  library(ggpubr)
  library(pheatmap)
  
  # Import BED files
  
  # Background: https://bioconductor.org/packages/devel/bioc/vignettes/BRGenomics/inst/doc/SpikeInNormalization.html
  
  files <- data.frame(
    uvb_rep1=c(
      'results/bed/IMR_200J_UVB_1_possort_markdup_spec1.bed.gz',
      'results/bed/IMR_200J_UVB_1_possort_markdup_spec2.bed.gz'
    ),
    uvb_rep2=c(
      'results/bed/IMR_200J_UVB_2_possort_markdup_spec1.bed.gz',
      'results/bed/IMR_200J_UVB_2_possort_markdup_spec2.bed.gz'
    ),
    uvc_rep1=c(
      'results/bed/IMR_100J_UVC_1_possort_markdup_spec1.bed.gz',
      'results/bed/IMR_100J_UVC_1_possort_markdup_spec2.bed.gz'
    ),
    uvc_rep2=c(
      'results/bed/IMR_100J_UVC_2_possort_markdup_spec1.bed.gz',
      'results/bed/IMR_100J_UVC_2_possort_markdup_spec2.bed.gz'
    )
  )
  
  
  
  # Import data and merge human with corresponding lambda counts
  bed <- apply(files,2,function(i) {
    mrg <- sapply(i,rtracklayer::import)
    mrg <- as(mrg,'GRangesList')
    suppressWarnings(mrg <- unlist(mrg))
    gc()
    unname(mrg)
  })
  bed <- GRangesList(bed,compress=TRUE)
  
  
  
  # Calculate normalization factors, 
  NF <- getSpikeInNFs(bed,si_names='chrL',method='SRPMC',batch_norm=TRUE,ctrl_pattern='uvb',field=NULL,ncores=1L)
  # Filter out spike reads
  bed.spike <- keepSeqlevels(bed,value='chrL',pruning.mode='fine')
  bed <- dropSeqlevels(bed,value='chrL',pruning.mode='fine')
  # Normalize counts
  tab <- lapply(bed,function(i) table(i$name))
  tab <- lapply(tab,'[',c('AA','AC','AG','AT','CA','CC','CG','CT','GA','GC','GG','GT','TA','TC','TG','TT'))
  tab <- do.call(rbind,tab)
  tab.norm <- as.data.frame(tab * NF)
  
  # saveRDS(bed,'method_bed.RDS')
  # saveRDS(bed.spike,'method_bed.spike.RDS')
  
  # Plot
  options(scipen=999)
  rn <- row.names(tab.norm)
  tab.norm$Library <- rn
  tab.norm$Condition <- sapply(strsplit(rn,'_'),'[',1)
  tab.norm$Replicate <- sapply(strsplit(rn,'_'),'[',2)
  # convert to long forat
  df <- pivot_longer(tab.norm,cols=where(is.numeric),names_to='Dinucleotide',values_to='SRPMC')
  # add summary statistics for error bars
  res <- df %>%
    group_by(Condition,Dinucleotide) %>%
    summarize(
      n=n(),
      mean=mean(SRPMC),
      sd=sd(SRPMC)
    ) %>%
    mutate(se=sd/sqrt(n))
  # do plotting
  pdf(file='figures/method.pdf',width=2,height=1.4)
  ggplot(res,aes(x=Dinucleotide,y=mean,fill=Condition)) + 
    geom_bar(stat='identity',position=position_dodge(width=0.9)) +
    geom_errorbar(aes(ymin=mean-se,ymax=mean+se),position=position_dodge(width=0.9),width=0.2) +
    #scale_fill_brewer(palette='Set1') +
    scale_fill_manual(values=c('#b2df8a','#a6cee3')) +
    guides(x=guide_axis(angle=90)) +
    theme(legend.position=c(0.5,1),
          legend.justification=c(1,1),
          panel.background=element_blank(),
          text=element_text(size=7))
  dev.off()
  
  
  
  # Comparison with Douki
  # https://pubs.acs.org/doi/abs/10.1021/bi0022543
  
  isolated_DNA_uvc <- data.frame(iso_uvc=c(2.970,1.823,0.573,0.069),Dinucleotide=c('TT','TC','CT','CC'))
  isolated_DNA_uvb <- data.frame(iso_uvb=c(1.023,0.694,0.289,0.094),Dinucleotide=c('TT','TC','CT','CC'))
  cellular_DNA_uvb <- data.frame(cel_uvb=c(3.147,1.286,0.577,0.279),Dinucleotide=c('TT','TC','CT','CC'))
  uvb <- res[res$Condition=='uvb',]
  uvb$TRAD_uvb <- uvb$mean
  uvc <- res[res$Condition=='uvc',]
  uvc$TRAD_uvc <- uvc$mean
  dat <- merge(merge(cellular_DNA_uvb,isolated_DNA_uvb,by='Dinucleotide'),isolated_DNA_uvc,by='Dinucleotide')
  dat <- merge(merge(dat,uvb,by='Dinucleotide',),uvc,by='Dinucleotide')

  #pdf(file='douki_compare.pdf',width=2.4,height=1.2)
  ggplot(dat,aes(x=cel_uvb,y=TRAD_uvb)) +
    geom_point(aes(shape=Dinucleotide)) +
    geom_line(linetype='dashed') +
    xlab('HPLC - MS/MS') +
    ylab('TRAD-Seq (SRPMC)') +
    ggtitle('Cellular DNA + UV-B') +
    stat_cor(method='pearson',digits=3) +
    theme(text=element_text(size=7))
  #dev.off()
  
  
  
  # Correlation heatmap
  # Normalize spike-in counts
  tab.spike <- lapply(bed.spike,function(i) table(i$name))
  tab.spike <- lapply(tab.spike,'[',c('AA','AC','AG','AT','CA','CC','CG','CT','GA','GC','GG','GT','TA','TC','TG','TT'))
  tab.spike <- do.call(rbind,tab.spike)
  tab.spike.norm <- as.data.frame(tab.spike * NF)
  rn.spike <- row.names(tab.spike.norm)
  tab.spike.norm$Library <- rn.spike
  tab.spike.norm$Condition <- sapply(strsplit(rn.spike,'_'),'[',1)
  tab.spike.norm$Replicate <- sapply(strsplit(rn.spike,'_'),'[',2)
  df.spike <- pivot_longer(tab.spike.norm,cols=where(is.numeric),names_to='Dinucleotide',values_to='SRPMC')
  # add summary statistics for error bars
  res.spike <- df.spike %>%
    group_by(Condition,Dinucleotide) %>%
    summarize(
      n=n(),
      mean=mean(SRPMC),
      sd=sd(SRPMC)
    ) %>%
    mutate(se=sd/sqrt(n))
  
  mat <- cbind(dat[,c('cel_uvb','iso_uvb','iso_uvc','TRAD_uvb','TRAD_uvc')],
               TRAD_iso_uvb=filter(res.spike,Condition=='uvb',Dinucleotide %in% c('TT','TC','CT','CC'))$mean,
               TRAD_iso_uvc=filter(res.spike,Condition=='uvc',Dinucleotide %in% c('TT','TC','CT','CC'))$mean)
  colnames(mat) <- c('HPLC_cellular_UV-B','HPLC_isolated_UV-B','HPLC_isolated_UV-C','TRAD_cellular_UV-B','TRAD_cellular_UV-C','TRAD_lambda_UV-B','TRAD_lambda_UV-C')
  corMat <- cor(mat,method='pearson')
  
  cols <- rev(brewer.pal(9,'Spectral'))
  pheatmap(corMat,display_numbers=TRUE,breaks=seq(0.8,1,length.out=10),color=cols)
  
  
  
  
  
  
  
  
  
  
