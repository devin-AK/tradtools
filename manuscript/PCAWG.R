
  setwd('~/Desktop/trad/PCAWG/manuscript/')
  
  library(R.utils)
  library(data.table)
  
  tcga <- fread('final_consensus_passonly.snv_mnv_indel.tcga.controlled.maf.gz')
  icgc <- fread('final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz')
  
  tmp <- table(tcga$Project_Code)
  
  barplot(tmp,las=2)
  
  tmp2 <- tmp[tmp$Donor_ID=='D046416',]