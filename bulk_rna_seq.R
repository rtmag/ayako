library(Rsubread)
options(scipen=999)

data<-featureCounts(c(
"/root/ayako/bam/Ctrl_CD41minus1_Aligned.sortedByCoord.out.bam",
"/root/ayako/bam/Ctrl_CD41minus2_Aligned.sortedByCoord.out.bam",
"/root/ayako/bam/Ctrl_CD41minus3_Aligned.sortedByCoord.out.bam",
"/root/ayako/bam/Ctrl_CD41plus1_Aligned.sortedByCoord.out.bam",
"/root/ayako/bam/Ctrl_CD41plus2_Aligned.sortedByCoord.out.bam",
"/root/ayako/bam/Ctrl_CD41plus3_Aligned.sortedByCoord.out.bam",
"/root/ayako/bam/Thpo_CD41minus1_Aligned.sortedByCoord.out.bam",
"/root/ayako/bam/Thpo_CD41minus2_Aligned.sortedByCoord.out.bam",
"/root/ayako/bam/Thpo_CD41minus3_Aligned.sortedByCoord.out.bam",
"/root/ayako/bam/Thpo_CD41plus1_Aligned.sortedByCoord.out.bam",
"/root/ayako/bam/Thpo_CD41plus2_Aligned.sortedByCoord.out.bam",
"/root/ayako/bam/Thpo_CD41plus3_Aligned.sortedByCoord.out.bam",
"/root/ayako/bam/Ctrl_neg_Aligned.sortedByCoord.out.bam"),
annot.ext="/root/ayako/ref/gencode.vM15.annotation.gtf",
isGTFAnnotationFile=TRUE,
minMQS=10,
strandSpecific=1,
isPairedEnd=TRUE,
autosort=TRUE,
nthreads=40,
GTF.attrType="gene_name"
)

dat=data[[1]]
colnames(dat)<-c("Ctrl_CD41minus1","Ctrl_CD41minus2","Ctrl_CD41minus3","Ctrl_CD41plus1","Ctrl_CD41plus2","Ctrl_CD41plus3",
                 "Thpo_CD41minus1","Thpo_CD41minus2","Thpo_CD41minus3","Thpo_CD41plus1","Thpo_CD41plus2","Thpo_CD41plus3",
                 "Ctrl_neg")
saveRDS(dat,"ayako_bulk_rna_counts.rds")
#################################################################################################
# DESEQ2
countData=readRDS("ayako_bulk_rna_counts.rds")
library(Rsubread)
options(scipen=999)
library(DESeq2)
#
design<-data.frame(group=c("Ctrl_CD41minus","Ctrl_CD41minus","Ctrl_CD41minus",
                           "Ctrl_CD41plus","Ctrl_CD41plus","Ctrl_CD41plus",
                           "Thpo_CD41minus","Thpo_CD41minus","Thpo_CD41minus",
                           "Thpo_CD41plus","Thpo_CD41plus","Thpo_CD41plus",
                           "Ctrl_neg") )


dLRT <- DESeqDataSetFromMatrix(countData = countData, colData = design, design = ~ group )
dLRT <- DESeq(dLRT, test="LRT", reduced=~1)
dLRT_vsd <- varianceStabilizingTransformation(dLRT)
vsd = assay(dLRT_vsd)

pdf("Diagnostic_pca.pdf")
plotPCA(dLRT_vsd,ntop=50000,intgroup=c("group"))
dev.off()

# rm2
countData=readRDS("ayako_bulk_rna_counts.rds")
library(Rsubread)
options(scipen=999)
library(DESeq2)
#
design<-data.frame(group=c("Ctrl_CD41minus","Ctrl_CD41minus",
                           "Ctrl_CD41plus","Ctrl_CD41plus","Ctrl_CD41plus",
                           "Thpo_CD41minus","Thpo_CD41minus","Thpo_CD41minus",
                           "Thpo_CD41plus","Thpo_CD41plus","Thpo_CD41plus"
                           ) )


dLRT <- DESeqDataSetFromMatrix(countData = countData[,c(1:2,4:12)], colData = design, design = ~ group )
dLRT <- DESeq(dLRT, test="LRT", reduced=~1)
dLRT_vsd <- varianceStabilizingTransformation(dLRT)
vsd = assay(dLRT_vsd)

pdf("Diagnostic_pca_rm2.pdf")
plotPCA(dLRT_vsd,ntop=50000,intgroup=c("group"))
dev.off()



# Boxplotter

boxploter=function(name_gene){
bdata= vsd[rownames(countData) %in% name_gene]
  
bdata= data.frame( names = c(rep("CD41-_Ctrl",3),
                     rep("CD41+_Ctrl",3),
                     rep("CD41-_Thpo",3),
                     rep("CD41+_Thpo",3),
                     "Ctrl_neg" ) ,
                    values= bdata)
  
bdata$names<-factor(bdata$names, levels=c("CD41+_Ctrl", "CD41+_Thpo", "CD41-_Ctrl","CD41-_Thpo","Ctrl_neg"))
  
boxplot( values ~ names, data = bdata,
        main = name_gene, col=c("#ffb3ba","#ffdfba","#baffc9","#bae1ff","#ffffba"),outline=F,las=2,cex.names=.5)
  stripchart(values ~ names, vertical = TRUE, data = bdata, 
    method = "jitter", add = TRUE, pch = 20, col = 'red')
} 


pdf("boxplot_cell_cycle_genes.pdf")
par(mfrow=c(3,3),cex.axis=.80)
boxploter("Ccna1")
boxploter("Ccnb1")
boxploter("Ccne1")

boxploter("Cdkn1a")
boxploter("Cdkn1b")
boxploter("Cdkn1c")

boxploter("Cdkn2d")
boxploter("Myc")
boxploter("Myb")
dev.off()
###################################


pdf("boxplot_major_HSC_genes.pdf")
par(mfrow=c(3,3),cex.axis=.80)
boxploter("Kit")
boxploter("Procr")
boxploter("Mapk14")

boxploter("Slamf1")
boxploter("Akt1")
boxploter("Cdh2")

boxploter("Atm")
boxploter("Atr")
boxploter("Egr1")
#
boxploter("Egr2")
boxploter("Cbfa2t3")
boxploter("Cd34")

boxploter("Fancc")
boxploter("Fanca")
boxploter("Fbxw7")

boxploter("Flt3")
boxploter("Gfi1")
boxploter("Gfi1b")
#
boxploter("Hif1a")
boxploter("Meis1")
boxploter("Hoxa9")

boxploter("Hoxb4")
boxploter("Itga4")
boxploter("Itga9")

boxploter("Tal1")
boxploter("Tek")
boxploter("Etv6")
#
boxploter("Bmi1")
boxploter("Hoxb5")
boxploter("Mecom")

boxploter("Adgrg1")
dev.off()

###################################

pdf("boxplot_major_MK_genes.pdf")
par(mfrow=c(3,3),cex.axis=.80)
boxploter("Mkl1")
boxploter("Mpl")
boxploter("Nfe2")

boxploter("Pf4")
boxploter("Thpo")
boxploter("Vwf")

boxploter("Zfpm1")
boxploter("Zfpm2")
boxploter("Itga2b")
#
boxploter("Cd9")
boxploter("Eng")
boxploter("Fli1")

boxploter("Runx1")
boxploter("Runx2")
boxploter("Runx3")

boxploter("Gata1")
boxploter("Gata2")
boxploter("Gata3")
#
boxploter("Flna")
boxploter("Pecam1")
boxploter("Pf4")
dev.off()
#
pdf("boxplot_TPO_signal_genes.pdf")
par(mfrow=c(3,3),cex.axis=.80)
boxploter("Jak2")
boxploter("Ptpn6")
boxploter("Ptpn11")

boxploter("Socs1")
boxploter("Socs2")
boxploter("Socs3")

boxploter("Stat1")
boxploter("Stat2")
boxploter("Stat3")
#
boxploter("Stat4")
boxploter("Stat5a")
boxploter("Stat5b")

boxploter("Cdc42")
boxploter("Gp1ba")
boxploter("Gp1bb")

boxploter("Gp5")
boxploter("Gp6")
boxploter("Gp9")
dev.off()
#

pdf("boxplot_mitochondrial_genes.pdf")
par(mfrow=c(3,3),cex.axis=.80)
boxploter("Atg12")
boxploter("Atg5")
boxploter("Atg7")

boxploter("Atp5a1")
boxploter("Cox5a")
boxploter("Cyc1")

boxploter("Cycs")
boxploter("Dnm1")
boxploter("Dnm1l")
#
boxploter("Flcn")
boxploter("Foxo3")
boxploter("Hif1a")

boxploter("Ldha")
boxploter("Map1lc3a")
boxploter("Mff")

boxploter("Mfn1")
boxploter("Mfn2")
boxploter("Mtor")
#
boxploter("Nfe2l2")
boxploter("Nrf1")
boxploter("Opa1")

boxploter("Park2")
boxploter("Pdk2")
boxploter("Pdk4")

boxploter("Pink1")
boxploter("Pparg")
boxploter("Ppargc1a")
#
boxploter("Sdha")
boxploter("Sdhb")
boxploter("Sdhc")

boxploter("Tfam")
boxploter("Tfe3")
boxploter("Ldb1")

boxploter("Mitf")
boxploter("Nfatc1")
boxploter("Pik3r1")
#
dev.off()

#################################################################################################

countData=readRDS("ayako_bulk_rna_counts.rds")

library(Rsubread)
options(scipen=999)
library(DESeq2)
library(gplots)
library(factoextra)
library(RColorBrewer)

countData = countData[,1:12]

design<-data.frame(group=c("Ctrl_CD41minus","Ctrl_CD41minus","Ctrl_CD41minus",
                           "Ctrl_CD41plus","Ctrl_CD41plus","Ctrl_CD41plus",
                           "Thpo_CD41minus","Thpo_CD41minus","Thpo_CD41minus",
                           "Thpo_CD41plus","Thpo_CD41plus","Thpo_CD41plus"
                           ) )


dLRT <- DESeqDataSetFromMatrix(countData = countData, colData = design, design = ~ group )
dLRT <- DESeq(dLRT, test="LRT", reduced=~1)
dLRT_vsd <- varianceStabilizingTransformation(dLRT)
vsd = assay(dLRT_vsd)
saveRDS(vsd,"vsd_bulk_ayako.rds")
dLRT_res = results(dLRT)
saveRDS(dLRT_res,"dLRT_res_bulk_ayako.rds")

vsd = readRDS("vsd_bulk_ayako.rds")
dLRT_res = readRDS("dLRT_res_bulk_ayako.rds")


postscript("anova.ps",height=10,width=10,horizontal=F)
sig_vsd = vsd[which(dLRT_res$padj<0.05),]
colnames(sig_vsd) <- c("CD41-_Ctrl","CD41-_Ctrl","CD41-_Ctrl",
                      "CD41+_Ctrl","CD41+_Ctrl","CD41+_Ctrl",
                      "CD41-_Thpo","CD41-_Thpo","CD41-_Thpo",
                      "CD41+_Thpo","CD41+_Thpo","CD41+_Thpo")
colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
heatmap.2(sig_vsd,col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="Genes",key.title="Gene expression",cexCol=.8)
dev.off()


postscript("anova2.ps",height=10,width=10,horizontal=F)
sig_vsd = vsd[which(dLRT_res$padj<0.05),]
colnames(sig_vsd) <- c("CD41-_Ctrl","CD41-_Ctrl","CD41-_Ctrl",
                      "CD41+_Ctrl","CD41+_Ctrl","CD41+_Ctrl",
                      "CD41-_Thpo","CD41-_Thpo","CD41-_Thpo",
                      "CD41+_Thpo","CD41+_Thpo","CD41+_Thpo")
colors <- colorRampPalette(c("blue","white","red"))(45)
heatmap.2(sig_vsd,col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="Genes",key.title="Gene expression",cexCol=.8)
dev.off()
#######################################################################################################################
options(bitmapType="cairo")
library(wordcloud)
countData=readRDS("ayako_bulk_rna_counts.rds")
vsd = readRDS("vsd_bulk_ayako.rds")

# 1) CD41+ UNTR VS CD41- UNTR
# DESEQ2
design<-data.frame(group=c("Ctrl_CD41minus","Ctrl_CD41minus","Ctrl_CD41minus",
                           "Ctrl_CD41plus","Ctrl_CD41plus","Ctrl_CD41plus") )
dds <- DESeqDataSetFromMatrix(countData = countData[,1:6], colData = design, design = ~ group )
dds <- DESeq(dds)
dds_res = results(dds,contrast=c("group","Ctrl_CD41plus","Ctrl_CD41minus"))

# Volcano

pdf("1_CD41+_untr_VS_CD41-_untr/volcano.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41+ Ctrl / CD41- Ctrl )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
       -log10(dds_res$padj)[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
      col="red",pch=20)
legend("topright", paste("CD41+ Ctrl:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41- Ctrl:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n") 

dev.off()

pdf("1_CD41+_untr_VS_CD41-_untr/volcano_labels.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41+ Ctrl / CD41- Ctrl )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
       -log10(dds_res$padj)[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
      col="blue",pch=20)
legend("topright", paste("CD41+ Ctrl:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41- Ctrl:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n") 

x=dds_res[which(dds_res$padj<0.05),]
x1=head(x[order(x$log2FoldChange),],15)
x2=tail(x[order(x$log2FoldChange),],15)
points(x1$log2FoldChange,-log10(x1$padj),pch=20,col="red")
nc=wordlayout(x1$log2FoldChange,-log10(x1$padj),rownames(x1),cex=1)
text(nc[,1],nc[,2],label=rownames(x1),cex=.7)
points(x2$log2FoldChange,-log10(x2$padj),pch=20,col="red")
nc=wordlayout(x2$log2FoldChange,-log10(x2$padj),rownames(x2),cex=1)
text(nc[,1],nc[,2],label=rownames(x2),cex=.7)
dev.off()

# Heatmap
postscript("1_CD41+_untr_VS_CD41-_untr/anova.ps",height=10,width=10,horizontal=F)
sig_vsd = vsd[which(abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05),1:6]
colnames(sig_vsd) <- c("CD41-_Ctrl","CD41-_Ctrl","CD41-_Ctrl",
                      "CD41+_Ctrl","CD41+_Ctrl","CD41+_Ctrl")

colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
heatmap.2(sig_vsd,col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="Genes",key.title="Gene expression",cexCol=.8)
dev.off()
# MAplot
pdf("1_CD41+_untr_VS_CD41-_untr/maplot.pdf")
plotMA(dds_res)
dev.off()
# RankedListFor GSEA
up_reg = dds_res[ which(dds_res$log2FoldChange>0),]
up_reg = up_reg[ !is.na(up_reg$padj),]
up_reg_log=-log(up_reg$padj)
names(up_reg_log) = rownames(up_reg)

dw_reg = dds_res[ which(dds_res$log2FoldChange<0),]
dw_reg = dw_reg[ !is.na(dw_reg$padj),]
dw_reg_log=log(dw_reg$padj)
names(dw_reg_log) = rownames(dw_reg)

rankedlist = cbind(sort(c(up_reg_log,dw_reg_log),decreasing=T) )
rankedlist = data.frame(ensid=rownames(rankedlist), log10FDR=rankedlist)
write.table(rankedlist,"1_CD41+_untr_VS_CD41-_untr/genes_ranked_table_FCFDR.rnk", sep="\t", quote=F,col.names=F,row.names=F)
# Significant Results ordered by log2FC
csv_table = dds_res[which(dds_res$padj<0.05 & abs(dds_res$log2FoldChange)>1),]
csv_table = csv_table[order(csv_table$log2FoldChange),]
write.csv(csv_table,"1_CD41+_untr_VS_CD41-_untr/differentially_expressed_genes.csv")

#######################################################################################################################
# 2) CD41+ TRTD VS CD41- TRTD
# DESEQ2

design<-data.frame(group=c("Thpo_CD41minus","Thpo_CD41minus","Thpo_CD41minus",
                           "Thpo_CD41plus","Thpo_CD41plus","Thpo_CD41plus"
                           ) )
                   
dds <- DESeqDataSetFromMatrix(countData = countData[,7:12], colData = design, design = ~ group )
dds <- DESeq(dds)
dds_res = results(dds,contrast=c("group","Thpo_CD41plus","Thpo_CD41minus"))

# Volcano

pdf("2_CD41+_tr_VS_CD41-_tr/volcano.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41+ Treated / CD41- Treated )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
       -log10(dds_res$padj)[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
      col="red",pch=20)
legend("topright", paste("CD41+ Treated:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41- Treated:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n") 

dev.off()

pdf("2_CD41+_tr_VS_CD41-_tr/volcano_labels.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41+ Treated / CD41- Treated )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
       -log10(dds_res$padj)[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
      col="blue",pch=20)
legend("topright", paste("CD41+ Treated:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41- Treated:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n") 

x=dds_res[which(dds_res$padj<0.05),]
x1=head(x[order(x$log2FoldChange),],15)
x2=tail(x[order(x$log2FoldChange),],15)
points(x1$log2FoldChange,-log10(x1$padj),pch=20,col="red")
nc=wordlayout(x1$log2FoldChange,-log10(x1$padj),rownames(x1),cex=1)
text(nc[,1],nc[,2],label=rownames(x1),cex=.7)
points(x2$log2FoldChange,-log10(x2$padj),pch=20,col="red")
nc=wordlayout(x2$log2FoldChange,-log10(x2$padj),rownames(x2),cex=1)
text(nc[,1],nc[,2],label=rownames(x2),cex=.7)
dev.off()

# Heatmap
postscript("2_CD41+_tr_VS_CD41-_tr/anova.ps",height=10,width=10,horizontal=F)
sig_vsd = vsd[which(abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05),7:12]
colnames(sig_vsd) <- c("CD41-_Thpo","CD41-_Thpo","CD41-_Thpo",
"CD41+_Thpo","CD41+_Thpo","CD41+_Thpo")

colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
heatmap.2(sig_vsd,col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="Genes",key.title="Gene expression",cexCol=.8)
dev.off()
# MAplot
pdf("2_CD41+_tr_VS_CD41-_tr/maplot.pdf")
plotMA(dds_res)
dev.off()
# RankedListFor GSEA
up_reg = dds_res[ which(dds_res$log2FoldChange>0),]
up_reg = up_reg[ !is.na(up_reg$padj),]
up_reg_log=-log(up_reg$padj)
names(up_reg_log) = rownames(up_reg)

dw_reg = dds_res[ which(dds_res$log2FoldChange<0),]
dw_reg = dw_reg[ !is.na(dw_reg$padj),]
dw_reg_log=log(dw_reg$padj)
names(dw_reg_log) = rownames(dw_reg)

rankedlist = cbind(sort(c(up_reg_log,dw_reg_log),decreasing=T) )
rankedlist = data.frame(ensid=rownames(rankedlist), log10FDR=rankedlist)
write.table(rankedlist,"2_CD41+_tr_VS_CD41-_tr/genes_ranked_table_FCFDR.rnk", sep="\t", quote=F,col.names=F,row.names=F)
# Significant Results ordered by log2FC
csv_table = dds_res[which(dds_res$padj<0.05 & abs(dds_res$log2FoldChange)>1),]
csv_table = csv_table[order(csv_table$log2FoldChange),]
write.csv(csv_table,"2_CD41+_tr_VS_CD41-_tr/differentially_expressed_genes.csv")


#######################################################################################################################
# 3) CD41+ UNTR VS CD41+ TRTD
# DESEQ2
design<-data.frame(group=c("Ctrl_CD41plus","Ctrl_CD41plus","Ctrl_CD41plus",
                           "Thpo_CD41plus","Thpo_CD41plus","Thpo_CD41plus"
                           ) )

dds <- DESeqDataSetFromMatrix(countData = countData[,c(4:6,10:12)], colData = design, design = ~ group )
dds <- DESeq(dds)
dds_res = results(dds,contrast=c("group","Ctrl_CD41plus","Thpo_CD41plus"))
# Volcano

pdf("3_CD41+_untr_VS_CD41+_tr/volcano.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41+ Untreated / CD41+ Treated )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
       -log10(dds_res$padj)[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
      col="red",pch=20)
legend("topright", paste("CD41+ Untreated:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41+ Treated:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n") 

dev.off()

pdf("3_CD41+_untr_VS_CD41+_tr/volcano_labels.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41+ Untreated / CD41+ Treated )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
       -log10(dds_res$padj)[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
      col="blue",pch=20)
legend("topright", paste("CD41+ Untreated:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41+ Treated:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n") 

x=dds_res[which(dds_res$padj<0.05),]
x1=head(x[order(x$log2FoldChange),],15)
x2=tail(x[order(x$log2FoldChange),],15)
points(x1$log2FoldChange,-log10(x1$padj),pch=20,col="red")
nc=wordlayout(x1$log2FoldChange,-log10(x1$padj),rownames(x1),cex=1)
text(nc[,1],nc[,2],label=rownames(x1),cex=.7)
points(x2$log2FoldChange,-log10(x2$padj),pch=20,col="red")
nc=wordlayout(x2$log2FoldChange,-log10(x2$padj),rownames(x2),cex=1)
text(nc[,1],nc[,2],label=rownames(x2),cex=.7)
dev.off()

# Heatmap
postscript("3_CD41+_untr_VS_CD41+_tr/anova.ps",height=10,width=10,horizontal=F)
sig_vsd = vsd[which(abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05),c(4:6,10:12)]
colnames(sig_vsd) <- c("CD41+_Ctrl","CD41+_Ctrl","CD41+_Ctrl",
"CD41+_Thpo","CD41+_Thpo","CD41+_Thpo")

colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
heatmap.2(sig_vsd,col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="Genes",key.title="Gene expression",cexCol=.8)
dev.off()
# MAplot
pdf("3_CD41+_untr_VS_CD41+_tr/maplot.pdf")
plotMA(dds_res)
dev.off()
# RankedListFor GSEA
up_reg = dds_res[ which(dds_res$log2FoldChange>0),]
up_reg = up_reg[ !is.na(up_reg$padj),]
up_reg_log=-log(up_reg$padj)
names(up_reg_log) = rownames(up_reg)

dw_reg = dds_res[ which(dds_res$log2FoldChange<0),]
dw_reg = dw_reg[ !is.na(dw_reg$padj),]
dw_reg_log=log(dw_reg$padj)
names(dw_reg_log) = rownames(dw_reg)

rankedlist = cbind(sort(c(up_reg_log,dw_reg_log),decreasing=T) )
rankedlist = data.frame(ensid=rownames(rankedlist), log10FDR=rankedlist)
write.table(rankedlist,"3_CD41+_untr_VS_CD41+_tr/genes_ranked_table_FCFDR.rnk", sep="\t", quote=F,col.names=F,row.names=F)
# Significant Results ordered by log2FC
csv_table = dds_res[which(dds_res$padj<0.05 & abs(dds_res$log2FoldChange)>1),]
csv_table = csv_table[order(csv_table$log2FoldChange),]
write.csv(csv_table,"3_CD41+_untr_VS_CD41+_tr/differentially_expressed_genes.csv")

#######################################################################################################################
# 4) CD41- UNTR VS CD41- TRTD
# DESEQ2
design<-data.frame(group=c("Ctrl_CD41minus","Ctrl_CD41minus","Ctrl_CD41minus",
                           "Thpo_CD41minus","Thpo_CD41minus","Thpo_CD41minus"
                           ) )

dds <- DESeqDataSetFromMatrix(countData = countData[,c(1:3,7:9)], colData = design, design = ~ group )
dds <- DESeq(dds)
dds_res = results(dds,contrast=c("group","Ctrl_CD41minus","Thpo_CD41minus"))
# Volcano

pdf("4_CD41-_untr_VS_CD41-_tr/volcano.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41- Untreated / CD41- Treated )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
       -log10(dds_res$padj)[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
      col="red",pch=20)
legend("topright", paste("CD41- Untreated:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41- Treated:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n") 

dev.off()

pdf("4_CD41-_untr_VS_CD41-_tr/volcano_labels.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41- Untreated / CD41- Treated )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
       -log10(dds_res$padj)[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
      col="blue",pch=20)
legend("topright", paste("CD41- Untreated:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41- Treated:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n") 

x=dds_res[which(dds_res$padj<0.05),]
x1=head(x[order(x$log2FoldChange),],15)
x2=tail(x[order(x$log2FoldChange),],15)
points(x1$log2FoldChange,-log10(x1$padj),pch=20,col="red")
nc=wordlayout(x1$log2FoldChange,-log10(x1$padj),rownames(x1),cex=1)
text(nc[,1],nc[,2],label=rownames(x1),cex=.7)
points(x2$log2FoldChange,-log10(x2$padj),pch=20,col="red")
nc=wordlayout(x2$log2FoldChange,-log10(x2$padj),rownames(x2),cex=1)
text(nc[,1],nc[,2],label=rownames(x2),cex=.7)
dev.off()

# Heatmap
postscript("4_CD41-_untr_VS_CD41-_tr/anova.ps",height=10,width=10,horizontal=F)
sig_vsd = vsd[which(abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05),c(1:3,7:9)]
colnames(sig_vsd) <- c("CD41-_Ctrl","CD41-_Ctrl","CD41-_Ctrl",
"CD41-_Thpo","CD41-_Thpo","CD41-_Thpo")

colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
heatmap.2(sig_vsd,col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="Genes",key.title="Gene expression",cexCol=.8)
dev.off()
# MAplot
pdf("4_CD41-_untr_VS_CD41-_tr/maplot.pdf")
plotMA(dds_res)
dev.off()
# RankedListFor GSEA
up_reg = dds_res[ which(dds_res$log2FoldChange>0),]
up_reg = up_reg[ !is.na(up_reg$padj),]
up_reg_log=-log(up_reg$padj)
names(up_reg_log) = rownames(up_reg)

dw_reg = dds_res[ which(dds_res$log2FoldChange<0),]
dw_reg = dw_reg[ !is.na(dw_reg$padj),]
dw_reg_log=log(dw_reg$padj)
names(dw_reg_log) = rownames(dw_reg)

rankedlist = cbind(sort(c(up_reg_log,dw_reg_log),decreasing=T) )
rankedlist = data.frame(ensid=rownames(rankedlist), log10FDR=rankedlist)
write.table(rankedlist,"4_CD41-_untr_VS_CD41-_tr/genes_ranked_table_FCFDR.rnk", sep="\t", quote=F,col.names=F,row.names=F)
# Significant Results ordered by log2FC
csv_table = dds_res[which(dds_res$padj<0.05 & abs(dds_res$log2FoldChange)>1),]
csv_table = csv_table[order(csv_table$log2FoldChange),]
write.csv(csv_table,"4_CD41-_untr_VS_CD41-_tr/differentially_expressed_genes.csv")


#######################################################################################################################
# 5) CD41+ VS CD41- Controlled by Treatment
# DESEQ2
design<-data.frame(CD41=c("CD41minus","CD41minus","CD41minus",
                           "CD41plus","CD41plus","CD41plus",
                           "CD41minus","CD41minus","CD41minus",
                           "CD41plus","CD41plus","CD41plus"),
                   Treatment=c("Ctrl","Ctrl","Ctrl",
                           "Ctrl","Ctrl","Ctrl",
                           "Thpo","Thpo","Thpo",
                           "Thpo","Thpo","Thpo") )

dds <- DESeqDataSetFromMatrix(countData = countData[,1:12], colData = design, 
                  design = ~ Treatment + CD41 )

dds <- DESeq(dds, test="LRT", 
           full= ~ Treatment + CD41, 
           reduced= ~ Treatment )
dds_res <- results(dds,contrast=c("CD41","CD41plus","CD41minus"))

# Volcano

pdf("5_CD41+_VS_CD41-_ControlledByTreatment/volcano.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41+ / CD41- )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
       -log10(dds_res$padj)[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
      col="red",pch=20)
legend("topright", paste("CD41+:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n")
legend("topleft", paste("CD41-:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n")

dev.off()

pdf("5_CD41+_VS_CD41-_ControlledByTreatment/volcano_labels.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41+ / CD41- )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
       -log10(dds_res$padj)[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
      col="blue",pch=20)
legend("topright", paste("CD41+:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41-:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n") 

x=dds_res[which(dds_res$padj<0.05),]
x1=head(x[order(x$log2FoldChange),],15)
x2=tail(x[order(x$log2FoldChange),],15)
points(x1$log2FoldChange,-log10(x1$padj),pch=20,col="red")
nc=wordlayout(x1$log2FoldChange,-log10(x1$padj),rownames(x1),cex=1)
text(nc[,1],nc[,2],label=rownames(x1),cex=.7)
points(x2$log2FoldChange,-log10(x2$padj),pch=20,col="red")
nc=wordlayout(x2$log2FoldChange,-log10(x2$padj),rownames(x2),cex=1)
text(nc[,1],nc[,2],label=rownames(x2),cex=.7)
dev.off()

# Heatmap
postscript("5_CD41+_VS_CD41-_ControlledByTreatment/anova.ps",height=10,width=10,horizontal=F)
sig_vsd = vsd[which(abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05),]
colnames(sig_vsd) <- c("CD41-_Ctrl","CD41-_Ctrl","CD41-_Ctrl",
"CD41+_Ctrl","CD41+_Ctrl","CD41+_Ctrl",
"CD41-_Thpo","CD41-_Thpo","CD41-_Thpo",
"CD41+_Thpo","CD41+_Thpo","CD41+_Thpo")

colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
heatmap.2(sig_vsd,col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="Genes",key.title="Gene expression",cexCol=.8)
dev.off()
# MAplot
pdf("5_CD41+_VS_CD41-_ControlledByTreatment/maplot.pdf")
plotMA(dds_res)
dev.off()
# RankedListFor GSEA
up_reg = dds_res[ which(dds_res$log2FoldChange>0),]
up_reg = up_reg[ !is.na(up_reg$padj),]
up_reg_log=-log(up_reg$padj)
names(up_reg_log) = rownames(up_reg)

dw_reg = dds_res[ which(dds_res$log2FoldChange<0),]
dw_reg = dw_reg[ !is.na(dw_reg$padj),]
dw_reg_log=log(dw_reg$padj)
names(dw_reg_log) = rownames(dw_reg)

rankedlist = cbind(sort(c(up_reg_log,dw_reg_log),decreasing=T) )
rankedlist = data.frame(ensid=rownames(rankedlist), log10FDR=rankedlist)
write.table(rankedlist,"5_CD41+_VS_CD41-_ControlledByTreatment/genes_ranked_table_FCFDR.rnk", sep="\t", quote=F,col.names=F,row.names=F)
# Significant Results ordered by log2FC
csv_table = dds_res[which(dds_res$padj<0.05 & abs(dds_res$log2FoldChange)>1),]
csv_table = csv_table[order(csv_table$log2FoldChange),]
write.csv(csv_table,"5_CD41+_VS_CD41-_ControlledByTreatment/differentially_expressed_genes.csv")


#######################################################################################################################
# 6) UNTR VS TRTD Controlled by CD41 status
# DESEQ2

design<-data.frame(CD41=c("CD41minus","CD41minus","CD41minus",
                           "CD41plus","CD41plus","CD41plus",
                           "CD41minus","CD41minus","CD41minus",
                           "CD41plus","CD41plus","CD41plus"),
                   Treatment=c("Ctrl","Ctrl","Ctrl",
                           "Ctrl","Ctrl","Ctrl",
                           "Thpo","Thpo","Thpo",
                           "Thpo","Thpo","Thpo") )

dds <- DESeqDataSetFromMatrix(countData = countData[,1:12], colData = design, 
                  design = ~ CD41 + Treatment )

dds <- DESeq(dds, test="LRT", 
           full= ~ CD41 + Treatment, 
           reduced= ~ CD41 )
dds_res <- results(dds,contrast=c("Treatment","Ctrl","Thpo"))

# Volcano

pdf("6_untr_VS_tr_ControlledByCD41Status/volcano.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( Ctrl / Thpo )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
       -log10(dds_res$padj)[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
      col="red",pch=20)
legend("topright", paste("Ctrl+:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n")
legend("topleft", paste("Thpo:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n")

dev.off()

pdf("6_untr_VS_tr_ControlledByCD41Status/volcano_labels.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( Ctrl / Thpo )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
       -log10(dds_res$padj)[abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05],
      col="blue",pch=20)
legend("topright", paste("Ctrl:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("Thpo:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n") 

x=dds_res[which(dds_res$padj<0.05),]
x1=head(x[order(x$log2FoldChange),],15)
x2=tail(x[order(x$log2FoldChange),],15)
points(x1$log2FoldChange,-log10(x1$padj),pch=20,col="red")
nc=wordlayout(x1$log2FoldChange,-log10(x1$padj),rownames(x1),cex=1)
text(nc[,1],nc[,2],label=rownames(x1),cex=.7)
points(x2$log2FoldChange,-log10(x2$padj),pch=20,col="red")
nc=wordlayout(x2$log2FoldChange,-log10(x2$padj),rownames(x2),cex=1)
text(nc[,1],nc[,2],label=rownames(x2),cex=.7)
dev.off()

# Heatmap
postscript("6_untr_VS_tr_ControlledByCD41Status/anova.ps",height=10,width=10,horizontal=F)
sig_vsd = vsd[which(abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05),]
colnames(sig_vsd) <- c("CD41-_Ctrl","CD41-_Ctrl","CD41-_Ctrl",
"CD41+_Ctrl","CD41+_Ctrl","CD41+_Ctrl",
"CD41-_Thpo","CD41-_Thpo","CD41-_Thpo",
"CD41+_Thpo","CD41+_Thpo","CD41+_Thpo")

colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
heatmap.2(sig_vsd,col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="Genes",key.title="Gene expression",cexCol=.8)
dev.off()
# MAplot
pdf("6_untr_VS_tr_ControlledByCD41Status/maplot.pdf")
plotMA(dds_res)
dev.off()
# RankedListFor GSEA
up_reg = dds_res[ which(dds_res$log2FoldChange>0),]
up_reg = up_reg[ !is.na(up_reg$padj),]
up_reg_log=-log(up_reg$padj)
names(up_reg_log) = rownames(up_reg)

dw_reg = dds_res[ which(dds_res$log2FoldChange<0),]
dw_reg = dw_reg[ !is.na(dw_reg$padj),]
dw_reg_log=log(dw_reg$padj)
names(dw_reg_log) = rownames(dw_reg)

rankedlist = cbind(sort(c(up_reg_log,dw_reg_log),decreasing=T) )
rankedlist = data.frame(ensid=rownames(rankedlist), log10FDR=rankedlist)
write.table(rankedlist,"6_untr_VS_tr_ControlledByCD41Status/genes_ranked_table_FCFDR.rnk", sep="\t", quote=F,col.names=F,row.names=F)
# Significant Results ordered by log2FC
csv_table = dds_res[which(dds_res$padj<0.05 & abs(dds_res$log2FoldChange)>1),]
csv_table = csv_table[order(csv_table$log2FoldChange),]
write.csv(csv_table,"6_untr_VS_tr_ControlledByCD41Status/differentially_expressed_genes.csv")







