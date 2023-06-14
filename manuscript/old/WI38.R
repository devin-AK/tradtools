### OLD
#se <- readRDS('~/Desktop/trad/counts/se_100Kb_100Kb.RDS')
### define batches
#batch3 <- c('1','2','1','2','1','2','1','2','2','3','1','2','1','2','1','2','1','2','1','2','1','2','1','2','1','2','1','3','2','3','2','3','2','3')
#batch2 <- ifelse(batch3=='1','A','B')
#batch4 <- paste(batch2,batch3,sep='_')
#colData(se) <- cbind(colData(se),batch2=as.factor(batch2),batch3=as.factor(batch3),batch4=as.factor(batch4))

# Estimate Size Factors from spike-in
#dds <- DESeqDataSet(se,design=~batch2+condition)
dds <- DESeqDataSet(se,design=~batch+condition)
dds$condition <- relevel(dds$condition,ref='WI38_200J_UVB') # set reference level
controls <- grepl('^chrL:',row.names(se))
dds <- estimateSizeFactors(dds,type='ratio',locfunc=median,controlGenes=controls) # can be ratio, poscounts, or iterate

# Differential susceptibility analysis
dds <- dds[!grepl('^chrL:',row.names(se))] # remove spike-in counts
dds <- DESeq(dds,parallel=TRUE)

Contrast <- c('condition','WI38RAF_200J_UVB','WI38_200J_UVB')
#Contrast <- c('condition','IMR_200J_UVB','WI38_200J_UVB')
res <- results(dds,contrast=Contrast)
sum(res$padj < 0.05,na.rm=TRUE)
sum(res$padj < 0.05,na.rm=TRUE) / nrow(res)
sum(res$padj < 0.05 & res$log2FoldChange>0,na.rm=T)
sum(res$padj < 0.05 & res$log2FoldChange<0,na.rm=T)

sum(res$padj < 0.00001,na.rm=TRUE)

resASH <- lfcShrink(dds, contrast=Contrast, type='ashr') # log fold-change shrinkage for visualizations
resASH$lambda_control <- grepl('^chrL:',row.names(resASH))
#resASH <- resASH[complete.cases(resASH),] # filtering
DSR_thresh <- c(0.05,0.01)
resASH$DSR <- ifelse(resASH$padj < DSR_thresh[1] & abs(resASH$log2FoldChange) > DSR_thresh[2],TRUE,FALSE)
resASH$minus_log10_pval <- -1*log10(resASH$pvalue)


# MA and Volcano plots (with spike-ins)
pm1 <- ggplot(as.data.frame(resASH),aes(x=baseMean,y=log2FoldChange)) +
  geom_point(aes(fill=lambda_control,color=DSR),pch=21) +
  scale_color_manual(values=c('gray','#FC9272')) +
  scale_fill_manual(values=c('gray','dodgerblue')) +
  # geom_hline(yintercept=c(-1,1),linetype='dashed') +
  scale_x_continuous(trans='log10',limits=c(0.01,15000)) +
  ylim(c(-4,4)) +
  labs(color='lambda control') +
  xlab('mean susceptibility') +
  ylab('log fold-change') +
  guides(color='none',fill='none') +
  theme_bw()
pv2 <- ggplot(as.data.frame(resASH),aes(x=log2FoldChange,y=minus_log10_pval)) +
  #geom_point(aes(color=DSR)) +
  geom_point(aes(fill=lambda_control,color=DSR),pch=21) +
  scale_color_manual(values=c('gray','#FC9272')) +
  scale_fill_manual(values=c('gray','dodgerblue')) +
  xlim(c(-4,4)) +
  ylim(c(0,20)) +
  ylab('-log10(pvalue)') +
  geom_hline(yintercept = -log10(DSR_thresh[1]), linetype='dashed') +
  #geom_vline(xintercept = c(-1*DSR_thresh[2],DSR_thresh[2]), linetype='dashed') + 
  theme_bw()
DSR_plots <- ggarrange(pm1,pv2,ncol=2,widths=c(1,1.5))
annotate_figure(DSR_plots,top=paste(Contrast[2],'vs',Contrast[3]))

# PCA
dds <- dds[,dds$condition %in% samples_to_plot]
dds$condition <- droplevels(dds$condition)
rld <- rlog(dds, blind=TRUE)
mat <- assay(rld)
mm <- model.matrix(~condition,colData(rld))
mat <- limma::removeBatchEffect(mat,batch=rld$batch,design=mm)
assay(rld) <- mat
ntop <- round(0.1*nrow(dds)) # top 10%
plotPCA(rld,ntop=ntop)

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
