cat CD41+_tr_fdr5_peaks.broadPeak CD41+_untr_fdr5_peaks.broadPeak CD41-_tr_fdr5_peaks.broadPeak| \
sort -k1,1 -k2,2n |bedtools merge -i - > CD41_merged_peaks.bed

#

bed_to_granges <- function(file){
   df <- read.table(file,
                    header=F,
                    stringsAsFactors=F)
 
   if(length(df) > 6){
      df <- df[,-c(7:length(df))]
   }
 
   if(length(df)<3){
      stop("File has less than 3 columns")
   }
 
   header <- c('chr','start','end','id','score','strand')
   names(df) <- header[1:length(names(df))]
 
   if('strand' %in% colnames(df)){
      df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
   }
 
   library("GenomicRanges")
 
   if(length(df)==3){
      gr <- with(df, GRanges(chr, IRanges(start, end)))
   } else if (length(df)==4){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
   } else if (length(df)==5){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
   } else if (length(df)==6){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
   }
   return(gr)
}

#
require(csaw)

param <- readParam(minq=10, pe='both')

regions=bed_to_granges("~/ayako/ayako_dejavu/peakcall/CD41_merged_peaks.bed")


counts <- regionCounts(bam.files, regions, param=param)
##
#
library(Rsubread)

x=read.table('~/ayako/ayako_dejavu/peakcall/CD41_merged_peaks.bed',sep="\t",stringsAsFactors=F)

ann = data.frame(GeneID=paste(x[,1],x[,2],x[,3],sep="_"),Chr=x[,1],Start=x[,2],End=x[,3],Strand='+')

bam.files <- c('/root/ayako/ayako_dejavu/bam/CD41+_untr_1_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41+_untr_2_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41+_untr_3_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41+_tr_1_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41+_tr_2_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41-_tr_1_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41-_tr_2_Aligned_rmdup.sortedByCoord.out.bam')



fc_SE <- featureCounts(bam.files,annot.ext=ann,isPairedEnd=TRUE,nthreads=20)

#
##

countData=fc_SE$counts

colnames(countData)=gsub('X.root.ayako.ayako_dejavu.bam.',"",colnames(countData))

colnames(countData)=gsub('_Aligned_rmdup.sortedByCoord.out.bam',"",colnames(countData))

##
#

countData=readRDS('atac_countdata.rds')


colnames(countData)=c("CD41_plus_untr_1","CD41_plus_untr_2","CD41_plus_untr_3","CD41_plus_tr_1","CD41_plus_tr_2","CD41_minus_tr_1","CD41_minus_tr_2")


require(DESeq2)

colData <- data.frame(group=colnames(countData))
dds <- DESeqDataSetFromMatrix(
       countData = countData,
       colData = colData,
       design = ~ group)

dLRT <- DESeq(dds, test="LRT", reduced=~1)
dLRT_vsd <- varianceStabilizingTransformation(dLRT)
dLRT_res <- results(dLRT)

pdf("Diagnostic_design_pca.pdf")
plotPCA(dLRT_vsd,ntop=30000,intgroup=c('group'))
dev.off()

####
# CD41+ untreated VS CD41+ treated

design<-data.frame(cells = c("CD41_plus_untr","CD41_plus_untr","CD41_plus_untr","CD41_plus_tr","CD41_plus_tr") )

dds <- DESeqDataSetFromMatrix(countData = countData[,c(1,2,3,4,5)], colData = design, design = ~ cells)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cells","CD41_plus_untr","CD41_plus_tr"))

library(graphics)
pdf("Volcano_CD41+_treated_vs_CD41+_untreated.pdf")
plot(res$log2FoldChange,-log10(res$padj),xlab=expression('Log'[2]*' Fold Change CD41+ UnTreated vs CD41+ Treated '),ylab=expression('-Log'[10]*' Q-values'),nrpoints=0)
dev.off()

ix=res$padj<0.05 & res$log2FoldChange<(-1)
ix[is.na(ix)]=FALSE
write.table(gsub("_","\t",rownames(res[ix,])),"CD41+_treated_over_CD41+_untreated.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)

ix=res$padj<0.05 & res$log2FoldChange>(1)
ix[is.na(ix)]=FALSE
write.table(gsub("_","\t",rownames(res[ix,])),"CD41+_untreated_over_CD41+_treated.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)
##

# CD41- treated VS CD41+ treated
design<-data.frame(cells = c("CD41_minus_tr","CD41_minus_tr","CD41_plus_tr","CD41_plus_tr") )

dds <- DESeqDataSetFromMatrix(countData = countData[,c(4,5,6,7)], colData = design, design = ~ cells)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cells","CD41_minus_tr","CD41_plus_tr"))

library(graphics)
pdf("Volcano_CD41-_treated_vs_CD41+_treated.pdf")
plot(res$log2FoldChange,-log10(res$padj),xlab=expression('Log'[2]*' Fold Change CD41- Treated vs CD41- Treated'),ylab=expression('-Log'[10]*' Q-values'),nrpoints=0)
dev.off()

ix=res$padj<0.05 & res$log2FoldChange<(-1)
ix[is.na(ix)]=FALSE
write.table(gsub("_","\t",rownames(res[ix,])),"CD41+_treated_over_CD41-_treated.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)

ix=res$padj<0.05 & res$log2FoldChange>(1)
ix[is.na(ix)]=FALSE
write.table(gsub("_","\t",rownames(res[ix,])),"CD41-_treated_over_CD41+_treated.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)

#############
# CD41+ untreated VS CD41- treated

design<-data.frame(cells = c("CD41_plus_untr","CD41_plus_untr","CD41_plus_untr","CD41_minus_tr","CD41_minus_tr") )

dds <- DESeqDataSetFromMatrix(countData = countData[,c(1,2,3,6,7)], colData = design, design = ~ cells)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cells","CD41_plus_untr","CD41_minus_tr"))

library(graphics)
pdf("Volcano_CD41+_untreated_vs_CD41-_treated.pdf")
plot(res$log2FoldChange,-log10(res$padj),xlab=expression('Log'[2]*' Fold Change CD41+ untreated vs CD41- Treated'),ylab=expression('-Log'[10]*' Q-values'),nrpoints=0)
dev.off()

ix=res$padj<0.05 & res$log2FoldChange<(-1)
ix[is.na(ix)]=FALSE
write.table(gsub("_","\t",rownames(res[ix,])),"CD41-_treated_over_CD41+_untreated.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)

ix=res$padj<0.05 & res$log2FoldChange>(1)
ix[is.na(ix)]=FALSE
write.table(gsub("_","\t",rownames(res[ix,])),"CD41+_untreated_over_CD41-_treated.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)
