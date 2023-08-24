
  setwd('~/Downloads')
  
  library(data.table)
  
  # c('CHEK2','CHEK1','ATR','ATM','TP53','MDM4','MDM2')
  
  # CHEK2
  chek2_fls <- dir('.',pattern='ENSG00000183765',full.names=T)
  
  # CHEK1
  chek1_fls <- dir('.',pattern='ENSG00000149554',full.names=T)
  
  # ATR
  atr_fls <- dir('.',pattern='ENSG00000175054',full.names=T)
  
  # ATM
  atm_fls <- dir('.',pattern='ENSG00000149311',full.names=T)
  
  # TP53
  tp53_fls <- dir('.',pattern='ENSG00000141510',full.names=T)
  
  # MDM4
  mdm4_fls <- dir('.',pattern='ENSG00000198625',full.names=T)
  
  # MDM2
  mdm2_fls <- dir('.',pattern='ENSG00000135679',full.names=T)
  
  
  
  # Process
  tmp <- function(x) {
    dat <- lapply(x,fread,data.table=FALSE)
    do.call('rbind',dat)
  }
  
  chek2 <- tmp(chek2_fls); chek2$Gene <- 'CHEK2'
  chek1 <- tmp(chek1_fls); chek1$Gene <- 'CHEK1'
  atr   <- tmp(atr_fls);   atr$Gene <- 'ATR'
  atm   <- tmp(atm_fls);   atm$Gene <- 'ATM'
  tp53  <- tmp(tp53_fls);  tp53$Gene <- 'TP53'
  mdm4  <- tmp(mdm4_fls);  mdm4$Gene <- 'MDM4'
  mdm2  <- tmp(mdm2_fls);  mdm2$Gene <- 'MDM2'
  
  x <- rbind(chek2,chek1,atr,atm,tp53,mdm4,mdm2)
  LET <- strsplit(x$`Genomic DNA Change`,split='\\d+')
  NUM <- strsplit(x$`Genomic DNA Change`,split='\\D+')
  
  x$Chromosome <- paste0(sapply(LET,'[',1),sapply(NUM,'[',2))
  x$Start <- sapply(NUM,'[',3)
  x$Mut <- sapply(LET,'[',3)
  x$Width <- sapply(strsplit(x$Mut,'>'),function(i)nchar(i[[1]]))
  
  write.table(x,'~/Desktop/trad/high_impact_p53_pathway_mutations.txt',quote=F,row.names=F,col.names=T,sep='\t')
  