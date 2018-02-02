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

colData <- data.frame(group=gsub("\\_\\d$","",colnames(countData),perl=TRUE) )
dds <- DESeqDataSetFromMatrix(
       countData = countData,
       colData = colData,
       design = ~ group)

dLRT <- DESeq(dds, test="LRT", reduced=~1)
dLRT_vsd <- varianceStabilizingTransformation(dLRT)
dLRT_res <- results(dLRT)

pdf("Diagnostic_design_pca.pdf")
plotPCA(dLRT_vsd,ntop=120000,intgroup=c('group'))
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

#################
#### intersect

intersectBed -a ~/ayako/ayako_dejavu/mm10_tss.bed -b CD41+_treated_over_CD41+_untreated.bed|cut -f4 > CD41+_treated_over_CD41+_untreated.tss
intersectBed -a ~/ayako/ayako_dejavu/mm10_tss.bed -b CD41+_untreated_over_CD41+_treated.bed|cut -f4 > CD41+_untreated_over_CD41+_treated.tss

intersectBed -a ~/ayako/ayako_dejavu/mm10_tss.bed -b CD41+_treated_over_CD41-_treated.bed|cut -f4 > CD41+_treated_over_CD41-_treated.tss
intersectBed -a ~/ayako/ayako_dejavu/mm10_tss.bed -b CD41-_treated_over_CD41+_treated.bed|cut -f4 > CD41-_treated_over_CD41+_treated.tss

intersectBed -a ~/ayako/ayako_dejavu/mm10_tss.bed -b CD41-_treated_over_CD41+_untreated.bed|cut -f4 > CD41-_treated_over_CD41+_untreated.tss
intersectBed -a ~/ayako/ayako_dejavu/mm10_tss.bed -b CD41+_untreated_over_CD41-_treated.bed|cut -f4 > CD41+_untreated_over_CD41-_treated.tss
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

legend('topright',legend=c("Open Chromatin in CD41+ treated","Open Chromatin in CD41+ untreated"),fill=c("#ffb3ba","#baffc9"),border=NA,bty = "n")
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

legend('topright',legend=c("Open Chromatin in CD41- treated","Open Chromatin in CD41+ treated"),fill=c("#bae1ff","#ffb3ba"),border=NA,bty = "n")
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

legend('topright',legend=c("Open Chromatin in CD41- treated","Open Chromatin in CD41+ untreated"),fill=c("#bae1ff","#baffc9"),border=NA,bty = "n")
dev.off()

################################################################################################
################################################
# ULTIMATE HEATMAP
source('https://raw.githubusercontent.com/rtmag/tumor-meth-pipe/master/heatmap3.R')
countData=readRDS('atac_countdata.rds')


colnames(countData)=c("CD41_plus_untr_1","CD41_plus_untr_2","CD41_plus_untr_3","CD41_plus_tr_1","CD41_plus_tr_2","CD41_minus_tr_1","CD41_minus_tr_2")


require(DESeq2)

colData <- data.frame(group=c("CD41_plus_untr","CD41_plus_untr","CD41_plus_untr","CD41_plus_tr","CD41_plus_tr","CD41_minus_tr","CD41_minus_tr"))
dds <- DESeqDataSetFromMatrix(
       countData = countData,
       colData = colData,
       design = ~ group)

dLRT <- DESeq(dds, test="LRT", reduced=~1)
dLRT_res <- results(dLRT)
dLRT_vsd <- varianceStabilizingTransformation(dLRT)
###FDR
#write.table(gsub("_","\t",rownames(dLRT_res[dLRT_res$padj<0.01 & !is.na(dLRT_res$padj),])),"LRT_FDR1.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)
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

design<-data.frame(cells = c("CD41_plus_untr","CD41_plus_untr","CD41_plus_untr","CD41_plus_tr","CD41_plus_tr") )

dds <- DESeqDataSetFromMatrix(countData = countData[,c(1,2,3,4,5)], colData = design, design = ~ cells)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cells","CD41_plus_untr","CD41_plus_tr"))
write.table(gsub("_","\t",rownames(res[res$padj<0.01 & !is.na(res$padj),])),"CD41+_trVSuntr_FDR1.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)

tss=read.table(pipe("intersectBed -a ../mm10_tss.bed -b CD41+_trVSuntr_FDR1.bed -wa -wb"),sep='\t',stringsAsFactors=F)
expr=read.table(pipe('grep -v "RNA-seq" ../GSE60101_1256271tableS2.txt'),sep="\t",header=T,stringsAsFactors=F)
vsd=assay(dLRT_vsd)

tss_vsd=vsd[match(paste(tss[,7],tss[,8],tss[,9],sep="_"),rownames(vsd)),1:5]

library(RColorBrewer)
colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(13)
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")

x=heatmap.3(tss_vsd,col=colors, hclustfun=hclustfunc, distfun=distfunc, 
            scale="row", trace="none",cexCol=1,KeyValueName="Expression",dendrogram="row")

match(tss[,4],expr[,1])


tss_vsd=vsd[rownames(vsd) %in% paste(tss[,7],tss[,8],tss[,9],sep="_"), 1:5]
x=heatmap.3(tss_vsd,col=colors,scale="row", trace="none",cexCol=1,KeyValueName="Expression",dendrogram="row")
##

# CD41- treated VS CD41+ treated
design<-data.frame(cells = c("CD41_minus_tr","CD41_minus_tr","CD41_plus_tr","CD41_plus_tr") )

dds <- DESeqDataSetFromMatrix(countData = countData[,c(4,5,6,7)], colData = design, design = ~ cells)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cells","CD41_minus_tr","CD41_plus_tr"))

write.table(gsub("_","\t",rownames(res[res$padj<0.05 & !is.na(res$padj),])),"CD41_+VS-_FDR5.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)

tss=read.table(pipe("intersectBed -a ../mm10_tss.bed -b CD41_+VS-_FDR5.bed -wa -wb"),sep='\t',stringsAsFactors=F)
expr=read.table(pipe('grep -v "RNA-seq" ../GSE60101_1256271tableS2.txt'),sep="\t",header=T,stringsAsFactors=F)
vsd=assay(dLRT_vsd)

tss_vsd=vsd[match(paste(tss[,7],tss[,8],tss[,9],sep="_"),rownames(vsd)),]

library(RColorBrewer)
colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(13)
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")

x=heatmap.3(tss_vsd[,4:7],col=colors, hclustfun=hclustfunc, distfun=distfunc, 
            scale="row", trace="none",cexCol=1,KeyValueName="Expression",dendrogram="row")

tss_vsd=vsd[rownames(vsd) %in% paste(tss[,7],tss[,8],tss[,9],sep="_"), 4:7]
x=heatmap.3(tss_vsd,col=colors,scale="row", trace="none",cexCol=1,KeyValueName="Expression",dendrogram="row")
##
########################################################################
##############################################################################
###########################################################################
ix = expr[,1] %in% tss[,4]
ex1=expr[ix, 3:18]
rownames(ex1)=make.names( expr[ix, 2], unique=T)




jx=tss[tss[,4] %in% expr[ix, 1],7:9]


####################################
#######################################

ix = expr[,1] %in% tss[,4]

ex1=expr[ix, 3:18]
rownames(ex1)=make.names( expr[ix, 2], unique=T)
ix
table(ix)
tss
table(ix)

ex1
 expr[ix, 1]
tss[,4] %in% expr[ix, 1]
table(tss[,4] %in% expr[ix, 1])
ix
table(tss[,4] %in% expr[ix, 1])
table(ix)
head(tss)
tss[tss[,4] %in% expr[ix, 1],7:9]
dim(tss[tss[,4] %in% expr[ix, 1],7:9])
ix
dim(tss[tss[,4] %in% expr[ix, 1],7:9])
table(ix)
dim(tss[tss[,4] %in% expr[ix, 1],7:9])
paste(tss[tss[,4] %in% expr[ix, 1],7:9])
paste(as.characters(tss[tss[,4] %in% expr[ix, 1],7:9]))
paste(as.character(tss[tss[,4] %in% expr[ix, 1],7:9]))
as.character(tss[tss[,4] %in% expr[ix, 1],7:9])
(tss[tss[,4] %in% expr[ix, 1],7:9])
paste(tss[tss[,4] %in% expr[ix, 1],7:9],sep="\t")
jx=tss[tss[,4] %in% expr[ix, 1],7:9]
jx
class(jx)
class(jx[,1])
class(jx[,2])
jx=tss[tss[,4] %in% expr[ix, 1],7:9]
tss=read.table(pipe("intersectBed -a ../mm10_tss.bed -b LRT_FDR5.bed -wa -wb"),sep='\t',stringsAsFactors=F)
expr=read.table(pipe('grep -v "RNA-seq" ../GSE60101_1256271tableS2.txt'),sep="\t",header=T,stringsAsFactors=F)
#
# LRT_FDR5
ix = expr[,1] %in% tss[,4]
ex1=expr[ix, 3:18]
rownames(ex1)=make.names( expr[ix, 2], unique=T)
jx=tss[tss[,4] %in% expr[ix, 1],7:9]
jx
class(jx[,2])
class(jx[,3])
class(jx[,1])
paste(jx)
paste(jx[,1])
paste(jx[,1],jx[,2])
paste(jx[,1],jx[,2],jx[,3])
head(paste(jx[,1],jx[,2],jx[,3]))
head(jx)
head(paste(jx[,1],jx[,2],jx[,3],sep="_"))
rown.name(dLRT_vsd) %in% (paste(jx[,1],jx[,2],jx[,3],sep="_")
rown.name(dLRT_vsd) %in% paste(jx[,1],jx[,2],jx[,3],sep="_")
row.name(dLRT_vsd) %in% paste(jx[,1],jx[,2],jx[,3],sep="_")
row.names(dLRT_vsd) %in% paste(jx[,1],jx[,2],jx[,3],sep="_")
table(row.names(dLRT_vsd) %in% paste(jx[,1],jx[,2],jx[,3],sep="_"))
tss
jx=tss[tss[,4] %in% expr[ix, 1],7:9]
jx
dim(jx)
dim(jx)
table(row.names(dLRT_vsd) %in% paste(jx[,1],jx[,2],jx[,3],sep="_"))
vsd[row.names(dLRT_vsd) %in% paste(jx[,1],jx[,2],jx[,3],sep="_"),]

