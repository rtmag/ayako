cat CD41+_tr_fdr5_peaks.broadPeak CD41+_untr_fdr5_peaks.broadPeak CD41-_tr_fdr5_peaks.broadPeak CD41-_untr_fdr5_peaks.broadPeak| \
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

ann = data.frame(GeneID=paste(x[,1],x[,2],x[,3],sep="_!_"),Chr=x[,1],Start=x[,2],End=x[,3],Strand='+')

bam.files <- c('/root/ayako/ayako_dejavu/bam/CD41+_untr_1_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41+_untr_2_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41+_untr_3_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41+_tr_1_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41+_tr_2_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41-_tr_1_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/bam/CD41-_tr_2_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/cd41-_untreated/bam/CD41-_untr_1_Aligned_rmdup.sortedByCoord.out.bam',
              '/root/ayako/ayako_dejavu/cd41-_untreated/bam/CD41-_untr_2_Aligned_rmdup.sortedByCoord.out.bam')




fc_SE <- featureCounts(bam.files,annot.ext=ann,isPairedEnd=TRUE,nthreads=20)

#
##

countData=fc_SE$counts

colnames(countData)=gsub('X.root.ayako.ayako_dejavu.bam.',"",colnames(countData))

colnames(countData)=gsub('_Aligned_rmdup.sortedByCoord.out.bam',"",colnames(countData))

saveRDS(countData,'atac_countdata.rds')
##
#

countData=readRDS('atac_countdata.rds')


colnames(countData)=c("CD41_plus_untr_1","CD41_plus_untr_2","CD41_plus_untr_3","CD41_plus_tr_1","CD41_plus_tr_2",
                      "CD41_minus_tr_1","CD41_minus_tr_2","CD41_minus_untr_1","CD41_minus_untr_2")


require(DESeq2)

colData <- data.frame(group=gsub("\\_\\d$","",colnames(countData),perl=TRUE) )
dds <- DESeqDataSetFromMatrix(
       countData = countData,
       colData = colData,
       design = ~ group)

dLRT <- DESeq(dds, test="LRT", reduced=~1)
dLRT_vsd <- varianceStabilizingTransformation(dLRT)
dLRT_res <- results(dLRT)
dLRT_res$padj[is.na(dLRT_res$padj)]=1


write.table(gsub("_!_","\t",rownames(dLRT_res[dLRT_res$padj<0.01,])),"ATAC-Seq_merged_LRT_FDR1.bed",
            quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(gsub("_!_","\t",rownames(dLRT_res[dLRT_res$padj<0.05,])),"ATAC-Seq_merged_LRT_FDR5.bed",
            quote=FALSE,col.names=FALSE,row.names=FALSE)


pdf("Diagnostic_design_pca.pdf")
plotPCA(dLRT_vsd,ntop=136500,intgroup=c('group'))
dev.off()

########################################
# CD41+ untreated VS CD41+ treated
design<-data.frame(cells = c("CD41_plus_untr","CD41_plus_untr","CD41_plus_untr","CD41_plus_tr","CD41_plus_tr") )

dds <- DESeqDataSetFromMatrix(countData = countData[,c(1,2,3,4,5)], colData = design, design = ~ cells)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cells","CD41_plus_untr","CD41_plus_tr"))

write.csv(res,"CD41+_untreated_vs_CD41+_treated.csv")

# CD41- untreated VS CD41- treated
design<-data.frame(cells = c("CD41_minus_untr","CD41_minus_untr","CD41_minus_tr","CD41_minus_tr") )

dds <- DESeqDataSetFromMatrix(countData = countData[,c(8,9,6,7)], colData = design, design = ~ cells)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cells","CD41_minus_untr","CD41_minus_tr"))

res_log=data.frame(res,log10=-log10(res$padj))

write.csv(res,"CD41-_untreated_vs_CD41-_treated.csv")
########################################
# CD41+ untreated VS CD41- untreated
design<-data.frame(cells = c("CD41_plus_untr","CD41_plus_untr","CD41_plus_untr","CD41_minus_untr","CD41_minus_untr") )

dds <- DESeqDataSetFromMatrix(countData = countData[,c(1,2,3,8,9)], colData = design, design = ~ cells)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cells","CD41_plus_untr","CD41_minus_untr"))

write.csv(res,"CD41+_untreated_vs_CD41-_untreated.csv")

# CD41+ treated VS CD41- treated
design<-data.frame(cells = c("CD41_plus_tr","CD41_plus_tr","CD41_minus_tr","CD41_minus_tr") )

dds <- DESeqDataSetFromMatrix(countData = countData[,c(4,5,6,7)], colData = design, design = ~ cells)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cells","CD41_plus_tr","CD41_minus_tr"))

write.csv(res,"CD41+_treated_vs_CD41-_treated.csv")

########################################
library(ggplot2)

res = read.csv("1_CD41+_untreated_vs_CD41+_treated.csv")

pdf("1_CD41+_untreated_vs_CD41+_treated.pdf")
plot(res$log2FoldChange,-log10(res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41+ Untreated / CD41+ Treated ) '),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.04))
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(res$log2FoldChange[abs(res$log2FoldChange)>1 & res$padj<0.05],
       -log10(res$padj)[abs(res$log2FoldChange)>1 & res$padj<0.05],
      col=alpha("#c0392b",.05))
legend("topright", paste("CD41+ untr:",length(which(res$log2FoldChange>1 & res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41+ tr:",length(which(res$log2FoldChange<(-1) & res$padj<0.05))), bty="n") 
dev.off()
#
res = read.csv("2_CD41-_untreated_vs_CD41-_treated.csv")

pdf("CD41-_untreated_vs_CD41-_treated.pdf")
plot(res$log2FoldChange,-log10(res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41- Untreated / CD41- Treated ) '),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.04))
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(res$log2FoldChange[abs(res$log2FoldChange)>1 & res$padj<0.05],
       -log10(res$padj)[abs(res$log2FoldChange)>1 & res$padj<0.05],
      col=alpha("#c0392b",.05))
legend("topright", paste("CD41- untr:",length(which(res$log2FoldChange>1 & res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41- tr:",length(which(res$log2FoldChange<(-1) & res$padj<0.05))), bty="n") 
dev.off()
#
res = read.csv("3_CD41+_untreated_vs_CD41-_untreated.csv")

pdf("3_CD41+_untreated_vs_CD41-_untreated.pdf")
plot(res$log2FoldChange,-log10(res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41+ Untreated / CD41- Untreated ) '),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.04))
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(res$log2FoldChange[abs(res$log2FoldChange)>1 & res$padj<0.05],
       -log10(res$padj)[abs(res$log2FoldChange)>1 & res$padj<0.05],
      col=alpha("#c0392b",.05))
legend("topright", paste("CD41+ untr:",length(which(res$log2FoldChange>1 & res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41- untr:",length(which(res$log2FoldChange<(-1) & res$padj<0.05))), bty="n") 
dev.off()
#
res = read.csv("4_CD41+_treated_vs_CD41-_treated.csv")

pdf("4_CD41+_treated_vs_CD41-_treated.pdf")
plot(res$log2FoldChange,-log10(res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41+ treated / CD41- Treated ) '),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.04))
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(res$log2FoldChange[abs(res$log2FoldChange)>1 & res$padj<0.05],
       -log10(res$padj)[abs(res$log2FoldChange)>1 & res$padj<0.05],
      col=alpha("#c0392b",.05))
legend("topright", paste("CD41+ tr:",length(which(res$log2FoldChange>1 & res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41- tr:",length(which(res$log2FoldChange<(-1) & res$padj<0.05))), bty="n")
dev.off()


##############################################################
res = read.csv("1_CD41+_untreated_vs_CD41+_treated.csv",row.name=1)
head(rownames(res[which(res$log2FoldChange>1 & res$padj<0.05),]))
CD41_plus_untr_1 = t(matrix(unlist(strsplit(rownames(res[which(res$log2FoldChange>1 & res$padj<0.05),]),"_!_")),nrow=3))
CD41_plus_tr_1 = t(matrix(unlist(strsplit(rownames(res[which(res$log2FoldChange<(-1) & res$padj<0.05),]),"_!_")),nrow=3))

write.table(CD41_plus_untr_1,"1_CD41+_untr_over_CD41+_tr.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(CD41_plus_tr_1,"1_CD41+_tr_over_CD41+_untr.bed",sep="\t",quote=F,row.names=F,col.names=F)

res = read.csv("2_CD41-_untreated_vs_CD41-_treated.csv",row.name=1)
CD41_minus_untr_2 = t(matrix(unlist(strsplit(rownames(res[which(res$log2FoldChange>1 & res$padj<0.05),]),"_!_")),nrow=3))
CD41_minus_tr_2 = t(matrix(unlist(strsplit(rownames(res[which(res$log2FoldChange<(-1) & res$padj<0.05),]),"_!_")),nrow=3))

write.table(CD41_minus_untr_2,"2_CD41-_untr_over_CD41-_tr.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(CD41_minus_tr_2,"2_CD41-_tr_over_CD41-_untr.bed",sep="\t",quote=F,row.names=F,col.names=F)

res = read.csv("3_CD41+_untreated_vs_CD41-_untreated.csv",row.name=1)
CD41_plus_untr_3 = t(matrix(unlist(strsplit(rownames(res[which(res$log2FoldChange>1 & res$padj<0.05),]),"_!_")),nrow=3))
CD41_minus_untr_3 = t(matrix(unlist(strsplit(rownames(res[which(res$log2FoldChange<(-1) & res$padj<0.05),]),"_!_")),nrow=3))

write.table(CD41_plus_untr_3,"3_CD41+_untr_over_CD41-_untr.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(CD41_minus_untr_3,"3_CD41-_untr_over_CD41+_untr.bed",sep="\t",quote=F,row.names=F,col.names=F)

res = read.csv("4_CD41+_treated_vs_CD41-_treated.csv",row.name=1)
CD41_plus_tr_4 = t(matrix(unlist(strsplit(rownames(res[which(res$log2FoldChange>1 & res$padj<0.05),]),"_!_")),nrow=3))
CD41_minus_tr_4 = t(matrix(unlist(strsplit(rownames(res[which(res$log2FoldChange<(-1) & res$padj<0.05),]),"_!_")),nrow=3))

write.table(CD41_plus_tr_4,"4_CD41+_tr_over_CD41-_tr.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(CD41_minus_tr_4,"4_CD41-_tr_over_CD41+_tr.bed",sep="\t",quote=F,row.names=F,col.names=F)

###############################################################
annotatePeaks.pl 1_CD41+_untr_over_CD41+_tr.bed mm10 > 1_CD41+_untr_over_CD41+_tr-mm10_annotation.txt
annotatePeaks.pl 1_CD41+_tr_over_CD41+_untr.bed mm10 > 1_CD41+_tr_over_CD41+_untr-mm10_annotation.txt

annotatePeaks.pl 2_CD41-_untr_over_CD41-_tr.bed mm10 > 2_CD41-_untr_over_CD41-_tr-mm10_annotation.txt
annotatePeaks.pl 2_CD41-_tr_over_CD41-_untr.bed mm10 > 2_CD41-_tr_over_CD41-_untr-mm10_annotation.txt

annotatePeaks.pl 3_CD41+_untr_over_CD41-_untr.bed mm10 > 3_CD41+_untr_over_CD41-_untr-mm10_annotation.txt
annotatePeaks.pl 3_CD41-_untr_over_CD41+_untr.bed mm10 > 3_CD41-_untr_over_CD41+_untr-mm10_annotation.txt

annotatePeaks.pl 4_CD41+_tr_over_CD41-_tr.bed mm10 > 4_CD41+_tr_over_CD41-_tr-mm10_annotation.txt
annotatePeaks.pl 4_CD41-_tr_over_CD41+_tr.bed mm10 > 4_CD41-_tr_over_CD41+_tr-mm10_annotation.txt

###############################################################

library(graphics)
library(ggplot2)

res_log=data.frame(res,log10=-log10(res$padj))

write.csv(res_log,"CD41+_treated_vs_CD41+_untreated.csv")

pdf("Volcano_CD41+_treated_vs_CD41+_untreated.pdf")
plot(res$log2FoldChange,-log10(res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41+ UnTreated / CD41+ Treated ) '),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("#f0650e",.5))
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
#f0650e
dev.off()


#ix=res$padj<0.05 & res$log2FoldChange<(-1)
#ix[is.na(ix)]=FALSE
#write.table(gsub("_","\t",rownames(res[ix,])),"results/CD41+_treated_over_CD41+_untreated.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)

#ix=res$padj<0.05 & res$log2FoldChange>(1)
#ix[is.na(ix)]=FALSE
#write.table(gsub("_","\t",rownames(res[ix,])),"results/CD41+_untreated_over_CD41+_treated.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)
##


# CD41- treated VS CD41+ treated
design<-data.frame(cells = c("CD41_minus_tr","CD41_minus_tr","CD41_plus_tr","CD41_plus_tr") )

dds <- DESeqDataSetFromMatrix(countData = countData[,c(4,5,6,7)], colData = design, design = ~ cells)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cells","CD41_minus_tr","CD41_plus_tr"))



#ix=res$padj<0.05 & res$log2FoldChange<(-1)
#ix[is.na(ix)]=FALSE
#write.table(gsub("_","\t",rownames(res[ix,])),"results/CD41+_treated_over_CD41-_treated.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)

#ix=res$padj<0.05 & res$log2FoldChange>(1)
#ix[is.na(ix)]=FALSE
#write.table(gsub("_","\t",rownames(res[ix,])),"results/CD41-_treated_over_CD41+_treated.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)

#############
# CD41+ untreated VS CD41- treated

design<-data.frame(cells = c("CD41_plus_untr","CD41_plus_untr","CD41_plus_untr","CD41_minus_tr","CD41_minus_tr") )

dds <- DESeqDataSetFromMatrix(countData = countData[,c(1,2,3,6,7)], colData = design, design = ~ cells)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cells","CD41_plus_untr","CD41_minus_tr"))

library(graphics)
pdf("Volcano_CD41+_untreated_vs_CD41-_treated.pdf")
plot(res$log2FoldChange,-log10(res$padj),xlab=expression('Log'[2]*' Fold Change CD41+ untreated vs CD41- Treated')
     ,ylab=expression('-Log'[10]*' Q-values'),nrpoints=0)
dev.off()

ix=res$padj<0.05 & res$log2FoldChange<(-1)
ix[is.na(ix)]=FALSE
write.table(gsub("_","\t",rownames(res[ix,])),"results/CD41-_treated_over_CD41+_untreated.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)

ix=res$padj<0.05 & res$log2FoldChange>(1)
ix[is.na(ix)]=FALSE
write.table(gsub("_","\t",rownames(res[ix,])),"results/CD41+_untreated_over_CD41-_treated.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)

#################
#### intersect

intersectBed -a ~/ayako/ayako_dejavu/mm10_tss.bed -b CD41+_treated_over_CD41+_untreated.bed|cut -f4 > \
CD41+_treated_over_CD41+_untreated.tss
intersectBed -a ~/ayako/ayako_dejavu/mm10_tss.bed -b CD41+_untreated_over_CD41+_treated.bed|cut -f4 > \
CD41+_untreated_over_CD41+_treated.tss

intersectBed -a ~/ayako/ayako_dejavu/mm10_tss.bed -b CD41+_treated_over_CD41-_treated.bed|cut -f4 > \
CD41+_treated_over_CD41-_treated.tss
intersectBed -a ~/ayako/ayako_dejavu/mm10_tss.bed -b CD41-_treated_over_CD41+_treated.bed|cut -f4 > \
CD41-_treated_over_CD41+_treated.tss

intersectBed -a ~/ayako/ayako_dejavu/mm10_tss.bed -b CD41-_treated_over_CD41+_untreated.bed|cut -f4 > \
CD41-_treated_over_CD41+_untreated.tss
intersectBed -a ~/ayako/ayako_dejavu/mm10_tss.bed -b CD41+_untreated_over_CD41-_treated.bed|cut -f4 > \
CD41+_untreated_over_CD41-_treated.tss
##
#
expr=read.table(pipe('grep -v "RNA-seq" ../GSE60101_1256271tableS2.txt'),sep="\t",header=T)
#
# CD41+ untreated VS CD41+ treated
x1=read.table('CD41+_treated_over_CD41+_untreated.tss')
x2=read.table('CD41+_untreated_over_CD41+_treated.tss')

ex1=expr[expr[,1] %in% x1[,1],3:18]
rownames(ex1)=make.names( expr[expr[,1] %in% x1[,1],2], unique=T)
ex2=expr[expr[,1] %in% x2[,1],3:18]
rownames(ex2)=make.names( expr[expr[,1] %in% x2[,1],2], unique=T)

track=c(rep(1,dim(ex1)[1]),rep(2,dim(ex2)[1]))

#ffb3ba red CD41+ treated
#baffc9 green CD41+ untreated
#bae1ff blue CD41- treated

colores=c("#ffb3ba","#baffc9")
rlab=rbind(colores[track])

 library(RColorBrewer)
colors <- colorRampPalette( (brewer.pal(9, "RdBu")) )(20)
cols=brewer.pal(3, "Set1")


# set the custom distance and clustering functions
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")

ex1.clust=hclustfunc(distfunc(ex1))
ex2.clust=hclustfunc(distfunc(ex2))

pdf("CD41+_untreated_VS_CD41+_treated.pdf")
x=heatmap.3(rbind(ex1[ex1.clust$order,],ex2[ex2.clust$order,]),col=colors, hclustfun=hclustfunc, distfun=distfunc, 
            scale="row", trace="none",cexCol=1,KeyValueName="Expression",Rowv=FALSE, RowSideColors=rlab,dendrogram="none")

legend('topright',legend=c("Open Chromatin in CD41+ treated","Open Chromatin in CD41+ untreated"),
       fill=c("#ffb3ba","#baffc9"),border=NA,bty = "n")
dev.off()
###

# CD41- treated VS CD41+ treated


x2=read.table('CD41+_treated_over_CD41-_treated.tss')
x1=read.table('CD41-_treated_over_CD41+_treated.tss')

ex1=expr[expr[,1] %in% x1[,1],3:18]
rownames(ex1)=make.names( expr[expr[,1] %in% x1[,1],2], unique=T)
ex2=expr[expr[,1] %in% x2[,1],3:18]
rownames(ex2)=make.names( expr[expr[,1] %in% x2[,1],2], unique=T)

track=c(rep(1,dim(ex1)[1]),rep(2,dim(ex2)[1]))

#ffb3ba red CD41+ treated
#baffc9 green CD41+ untreated
#bae1ff blue CD41- treated


colores=c("#bae1ff","#ffb3ba")
rlab=rbind(colores[track])

 library(RColorBrewer)
colors <- colorRampPalette( (brewer.pal(9, "RdBu")) )(20)


# set the custom distance and clustering functions
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")

ex1.clust=hclustfunc(distfunc(ex1))
ex2.clust=hclustfunc(distfunc(ex2))

pdf("CD41-_treated_VS_CD41+_treated.pdf")
x=heatmap.3(rbind(ex1,ex2[ex2.clust$order,]),col=colors, hclustfun=hclustfunc, distfun=distfunc, 
            scale="row", trace="none",cexCol=1,KeyValueName="Expression",Rowv=FALSE, RowSideColors=rlab,dendrogram="none")

legend('topright',legend=c("Open Chromatin in CD41- treated","Open Chromatin in CD41+ treated"),
       fill=c("#bae1ff","#ffb3ba"),border=NA,bty = "n")
dev.off()
###
# CD41+ untreated VS CD41- treated

x2=read.table('CD41+_untreated_over_CD41-_treated.tss')
x1=read.table('CD41-_treated_over_CD41+_untreated.tss')

ex1=expr[expr[,1] %in% x1[,1],3:18]
rownames(ex1)=make.names( expr[expr[,1] %in% x1[,1],2], unique=T)
ex2=expr[expr[,1] %in% x2[,1],3:18]
rownames(ex2)=make.names( expr[expr[,1] %in% x2[,1],2], unique=T)

track=c(rep(1,dim(ex1)[1]),rep(2,dim(ex2)[1]))

#ffb3ba red CD41+ treated
#baffc9 green CD41+ untreated
#bae1ff blue CD41- treated


colores=c("#bae1ff","#baffc9")
rlab=rbind(colores[track])

 library(RColorBrewer)
colors <- colorRampPalette( (brewer.pal(9, "RdBu")) )(20)


# set the custom distance and clustering functions
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")

ex1.clust=hclustfunc(distfunc(ex1))
ex2.clust=hclustfunc(distfunc(ex2))

pdf("CD41-_treated_VS_CD41+_untreated.pdf")
x=heatmap.3(rbind(ex1[ex1.clust$order,],ex2[ex2.clust$order,]),col=colors, hclustfun=hclustfunc, distfun=distfunc, 
            scale="row", trace="none",cexCol=1,KeyValueName="Expression",Rowv=FALSE, RowSideColors=rlab,dendrogram="none")

legend('topright',legend=c("Open Chromatin in CD41- treated","Open Chromatin in CD41+ untreated"),
       fill=c("#bae1ff","#baffc9"),border=NA,bty = "n")
dev.off()

################################################################################################
################################################
# ULTIMATE HEATMAP
source('https://raw.githubusercontent.com/rtmag/tumor-meth-pipe/master/heatmap3.R')
countData=readRDS('atac_countdata.rds')


colnames(countData)=c("CD41_plus_untr_1","CD41_plus_untr_2","CD41_plus_untr_3","CD41_plus_tr_1","CD41_plus_tr_2",
                      "CD41_minus_tr_1","CD41_minus_tr_2")


require(DESeq2)

colData <- data.frame(group=c("CD41_plus_untr","CD41_plus_untr","CD41_plus_untr","CD41_plus_tr","CD41_plus_tr",
                              "CD41_minus_tr","CD41_minus_tr"))
dds <- DESeqDataSetFromMatrix(
       countData = countData,
       colData = colData,
       design = ~ group)

dLRT <- DESeq(dds, test="LRT", reduced=~1)
dLRT_res <- results(dLRT)
dLRT_vsd <- varianceStabilizingTransformation(dLRT)
saveRDS(dLRT_res,"dLRT_res.rds")
saveRDS(dLRT_vsd,"dLRT_vsd.rds")

dLRT_res=readRDS('dLRT_res.rds')
dLRT_vsd=readRDS('dLRT_vsd.rds')
###FDR
#write.table(gsub("_","\t",rownames(dLRT_res[dLRT_res$padj<0.01 & !is.na(dLRT_res$padj),])),"LRT_FDR1.bed",
#quote=FALSE,col.names=FALSE,row.names=FALSE)
#
#tss=read.table(pipe("intersectBed -a ../mm10_tss.bed -b LRT_FDR5.bed -wa -wb"),sep='\t',stringsAsFactors=F)
tss=read.table(pipe("intersectBed -a ../mm10_tss.bed -b LRT_FDR1.bed -wa -wb"),sep='\t',stringsAsFactors=F)
expr=read.table(pipe('grep -v "RNA-seq" ../GSE60101_1256271tableS2.txt'),sep="\t",header=T,stringsAsFactors=F)
vsd=assay(dLRT_vsd)
#
tss_vsd=vsd[match(paste(tss[,7],tss[,8],tss[,9],sep="_"),rownames(vsd)),]

library(RColorBrewer)
colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(13)
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")

x=heatmap.3(tss_vsd,col=colors, hclustfun=hclustfunc, distfun=distfunc, 
            scale="row", trace="none",cexCol=1,KeyValueName="Expression",dendrogram="row")

##
tss_vsd=vsd[rownames(vsd) %in% paste(tss[,7],tss[,8],tss[,9],sep="_"), 1:5]
x=heatmap.3(tss_vsd[,1:5],col=colors,scale="row", trace="none",cexCol=1,KeyValueName="Expression",dendrogram="row")
##
colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(5)
tss_vsd=vsd[match(paste(tss[,7],tss[,8],tss[,9],sep="_"),rownames(vsd)),]
x=heatmap.3(tss_vsd,col=colors,scale="row", trace="none",cexCol=1,KeyValueName="Expression",dendrogram="row")

# CD41+ untreated VS CD41+ treated
dLRT_vsd=readRDS('dLRT_vsd.rds')

design<-data.frame(cells = c("CD41_plus_untr","CD41_plus_untr","CD41_plus_untr","CD41_plus_tr","CD41_plus_tr") )

dds <- DESeqDataSetFromMatrix(countData = countData[,c(1,2,3,4,5)], colData = design, design = ~ cells)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cells","CD41_plus_untr","CD41_plus_tr"))
write.table(gsub("_","\t",rownames(res[res$padj<0.01 & !is.na(res$padj),])),"CD41+_trVSuntr_FDR1.bed",
            quote=FALSE,col.names=FALSE,row.names=FALSE)

tss=read.table(pipe("intersectBed -a ../mm10_tss.bed -b CD41+_trVSuntr_FDR1.bed -wa -wb"),sep='\t',stringsAsFactors=F)
expr=read.table(pipe('grep -v "RNA-seq" ../GSE60101_1256271tableS2.txt'),sep="\t",header=T,stringsAsFactors=F)
vsd=assay(dLRT_vsd)

tss_vsd=vsd[match(paste(tss[,7],tss[,8],tss[,9],sep="_"),rownames(vsd)),1:5]

library(RColorBrewer)
colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(13)
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")
#distfunc <- function(x) as.dist(1-(cor(x, method="spearman")))

x=heatmap.3(tss_vsd,col=colors, hclustfun=hclustfunc, distfun=distfunc,
            scale="row", trace="none",cexRow=0.5,cexCol=.7,KeyValueName="Expression",dendrogram="row")
rowORDER=x$rowInd

tss_exp=expr[match(tss[,4],expr[,1]),]
tss_exp[is.na(tss_exp)] <- 0
rownames(tss_exp) = make.names( tss_exp[,2],unique=T)
tss_exp=(tss_exp[,3:18])
tss_exp=as.matrix(tss_exp)

colors <- colorRampPalette( (brewer.pal(9, "RdBu")) )(20)

heatmap.3(tss_exp[x$rowInd,],col=colors, hclustfun=hclustfunc, distfun=distfunc, 
            scale="row", trace="none",cexCol=1,KeyValueName="Expression",Rowv=FALSE,dendrogram="none")

##

# CD41- treated VS CD41+ treated
design<-data.frame(cells = c("CD41_minus_tr","CD41_minus_tr","CD41_plus_tr","CD41_plus_tr") )

dds <- DESeqDataSetFromMatrix(countData = countData[,c(4,5,6,7)], colData = design, design = ~ cells)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cells","CD41_minus_tr","CD41_plus_tr"))

write.table(gsub("_","\t",rownames(res[res$padj<0.05 & !is.na(res$padj),])),"CD41_+VS-_FDR5.bed",
            quote=FALSE,col.names=FALSE,row.names=FALSE)

tss=read.table(pipe("intersectBed -a ../mm10_tss.bed -b CD41_+VS-_FDR5.bed -wa -wb"),sep='\t',stringsAsFactors=F)
expr=read.table(pipe('grep -v "RNA-seq" ../GSE60101_1256271tableS2.txt'),sep="\t",header=T,stringsAsFactors=F)
vsd=assay(dLRT_vsd)

tss_vsd=vsd[match(paste(tss[,7],tss[,8],tss[,9],sep="_"),rownames(vsd)),4:7]

library(RColorBrewer)
colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(13)
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")

x=heatmap.3(tss_vsd,col=colors, hclustfun=hclustfunc, distfun=distfunc,
            scale="row", trace="none",cexRow=0.5,cexCol=.7,KeyValueName="Expression",dendrogram="row")
rowORDER=x$rowInd


tss_exp=expr[match(tss[,4],expr[,1]),]
tss_exp[is.na(tss_exp)] <- 0
rownames(tss_exp) = make.names( tss_exp[,2],unique=T)
tss_exp=(tss_exp[,3:18])
tss_exp=as.matrix(tss_exp)

colors <- colorRampPalette( (brewer.pal(9, "RdBu")) )(20)

#x=heatmap.3(tss_exp[x$rowInd,],col=colors, hclustfun=hclustfunc, distfun=distfunc, 
#            scale="row", trace="none",cexCol=1,KeyValueName="Expression",Rowv=FALSE,dendrogram="none")

##
########################################################################
##############################################################################
###########################################################################
library(gridGraphics)
library(grid)

library(gridGraphics)
grab_grob <- function(){
  grid.echo()
  grid.grab()
}

heatmap.3(tss_vsd,col=colors, hclustfun=hclustfunc, distfun=distfunc,
            scale="row", trace="none",cexRow=0.5,cexCol=.7,KeyValueName="Expression",dendrogram="row")
g1 <- grab_grob()

heatmap.3(tss_exp[rowORDER,],col=colors, hclustfun=hclustfunc, distfun=distfunc, 
            scale="col", trace="none",cexCol=1,KeyValueName="Expression",Rowv=FALSE,dendrogram="none")

g2 <- grab_grob()

pdf("CD41+_trVSuntr2.pdf")

grid.newpage()

# library(gridExtra)
# grid.arrange(g,g, ncol=2, clip=TRUE)

lay <- grid.layout(nrow = 1, ncol=2)
pushViewport(viewport(layout = lay))

grid.draw(editGrob(g1, vp=viewport(layout.pos.row = 1, 
                                  layout.pos.col = 1, clip=TRUE)))


grid.draw(editGrob(g2, vp=viewport(layout.pos.row = 1, 
                                  layout.pos.col = 2, clip=TRUE)))
upViewport(1)
dev.off()

########################################################################
##############################################################################
###########################################################################
heatmap.3(tss_vsd,col=colors, hclustfun=hclustfunc, distfun=distfunc,
            scale="row", trace="none",cexRow=0.5,cexCol=.7,KeyValueName="Expression",dendrogram="row")
g1 <- grab_grob()

x=heatmap.3(tss_exp[rowORDER,],col=colors, hclustfun=hclustfunc, distfun=distfunc, 
            scale="col", trace="none",cexCol=1,KeyValueName="Expression",Rowv=FALSE,dendrogram="none")

g2 <- grab_grob()

pdf("CD41+_tr_VS_CD41-_tr.pdf")

grid.newpage()

# library(gridExtra)
# grid.arrange(g,g, ncol=2, clip=TRUE)

lay <- grid.layout(nrow = 1, ncol=2)
pushViewport(viewport(layout = lay))

grid.draw(editGrob(g1, vp=viewport(layout.pos.row = 1, 
                                  layout.pos.col = 1, clip=TRUE)))


grid.draw(editGrob(g2, vp=viewport(layout.pos.row = 1, 
                                  layout.pos.col = 2, clip=TRUE)))
upViewport(1)
dev.off()

#############

intersectBed -a ../../mm10_tss.bed -b CD41+_treated_over_CD41+_untreated.bed -wa -wb| \
cut -f4| grep -f - ../../GSE60101_1256271tableS2.txt|cat header.txt - > CD41+_treated_over_CD41+_untreated_genes.txt
intersectBed -a ../../mm10_tss.bed -b CD41+_treated_over_CD41-_treated.bed -wa -wb| \
cut -f4| grep -f - ../../GSE60101_1256271tableS2.txt|cat header.txt - > CD41+_treated_over_CD41-_treated_genes.txt
intersectBed -a ../../mm10_tss.bed -b CD41+_untreated_over_CD41+_treated.bed -wa -wb| \
cut -f4| grep -f - ../../GSE60101_1256271tableS2.txt|cat header.txt - > CD41+_untreated_over_CD41+_treated_genes.txt
intersectBed -a ../../mm10_tss.bed -b CD41+_untreated_over_CD41-_treated.bed -wa -wb| \
cut -f4| grep -f - ../../GSE60101_1256271tableS2.txt|cat header.txt - > CD41+_untreated_over_CD41-_treated_genes.txt

intersectBed -a ../../mm10_tss.bed -b CD41-_treated_over_CD41+_treated.bed -wa -wb| \
cut -f4| grep -f - ../../GSE60101_1256271tableS2.txt|cat header.txt - > CD41-_treated_over_CD41+_treated_genes.txt
intersectBed -a ../../mm10_tss.bed -b CD41-_treated_over_CD41+_untreated.bed -wa -wb| \
cut -f4| grep -f - ../../GSE60101_1256271tableS2.txt|cat header.txt - > CD41-_treated_over_CD41+_untreated_genes.txt









