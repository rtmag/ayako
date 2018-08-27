
#######################################################################################################################
options(bitmapType="cairo")
library(wordcloud)
require(DESeq2)
library(gplots)
library(factoextra)
library(RColorBrewer)

countData=readRDS("ayako_bulk_rna_counts.rds")
vsd = readRDS("vsd_bulk_ayako.rds")

HR = c("Atm",
"Brca1",
"Bard1",
"Rbbp8",
"Rad50",
"Nbn",
"Rad51",
"Rad52",
"Exo1",
"Brca2",
"Palb2",
"Mre11a",
"Spo11",
"Mdc1",
"Rnf8",
"Xrcc2")

NHEJ = c("Trp53bp1",
         "Xrcc5",
"Xrcc6",
"Xrcc4",
"Nhej1",
"Ier3",
"Lig4",
"Prkdc",
"Wrn",
"Dclre1c")

# 1) CD41+ UNTR VS CD41- UNTR
# DESEQ2
design<-data.frame(group=c("Ctrl_CD41minus","Ctrl_CD41minus","Ctrl_CD41minus",
                           "Ctrl_CD41plus","Ctrl_CD41plus","Ctrl_CD41plus") )
dds <- DESeqDataSetFromMatrix(countData = countData[,1:6], colData = design, design = ~ group )
dds <- DESeq(dds)
dds_res = results(dds,contrast=c("group","Ctrl_CD41plus","Ctrl_CD41minus"))

# Volcano

pdf("1_CD41+_untr_VS_CD41-_untr_volcano_labels.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41+ Ctrl / CD41- Ctrl )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")

points(dds_res$log2FoldChange[rownames(dds_res) %in% HR],
       -log10(dds_res$padj)[rownames(dds_res) %in% HR],
      col="red",pch=20)

points(dds_res$log2FoldChange[rownames(dds_res) %in% NHEJ],
       -log10(dds_res$padj)[rownames(dds_res) %in% NHEJ],
      col="blue",pch=20)

legend("topright", paste("CD41+ Ctrl:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41- Ctrl:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n") 

x=dds_res
x1=x[which((rownames(dds_res) %in% NHEJ) & abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05), ]
x2=x[which((rownames(dds_res) %in% HR) & abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05), ]
nc=wordlayout(x1$log2FoldChange,-log10(x1$padj),rownames(x1),cex=1)
text(nc[,1],nc[,2],label=rownames(x1),cex=.7)
nc=wordlayout(x2$log2FoldChange,-log10(x2$padj),rownames(x2),cex=1)
text(nc[,1],nc[,2],label=rownames(x2),cex=.7)
legend("left", c("HR","NHEJ"),fill=c("red","blue"), bty="n") 

dev.off()


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

pdf("2_CD41+_tr_VS_CD41-_tr_volcano_labels.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41+ Treated / CD41- Treated )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[rownames(dds_res) %in% HR],
       -log10(dds_res$padj)[rownames(dds_res) %in% HR],
      col="red",pch=20)

points(dds_res$log2FoldChange[rownames(dds_res) %in% NHEJ],
       -log10(dds_res$padj)[rownames(dds_res) %in% NHEJ],
      col="blue",pch=20)
legend("topright", paste("CD41+ Treated:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41- Treated:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n") 

x=dds_res
x1=x[which((rownames(dds_res) %in% NHEJ) & abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05), ]
x2=x[which((rownames(dds_res) %in% HR) & abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05), ]
nc=wordlayout(x1$log2FoldChange,-log10(x1$padj),rownames(x1),cex=1)
text(nc[,1],nc[,2],label=rownames(x1),cex=.7)
nc=wordlayout(x2$log2FoldChange,-log10(x2$padj),rownames(x2),cex=1)
text(nc[,1],nc[,2],label=rownames(x2),cex=.7)
legend("left", c("HR","NHEJ"),fill=c("red","blue"), bty="n") 
dev.off()


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

pdf("3_CD41+_untr_VS_CD41+_tr/volcano_labels.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41+ Untreated / CD41+ Treated )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[rownames(dds_res) %in% HR],
       -log10(dds_res$padj)[rownames(dds_res) %in% HR],
      col="red",pch=20)

points(dds_res$log2FoldChange[rownames(dds_res) %in% NHEJ],
       -log10(dds_res$padj)[rownames(dds_res) %in% NHEJ],
      col="blue",pch=20)
legend("topright", paste("CD41+ Untreated:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41+ Treated:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n") 

x=dds_res
x1=x[which((rownames(dds_res) %in% NHEJ) & abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05), ]
x2=x[which((rownames(dds_res) %in% HR) & abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05), ]
nc=wordlayout(x1$log2FoldChange,-log10(x1$padj),rownames(x1),cex=1)
text(nc[,1],nc[,2],label=rownames(x1),cex=.7)
nc=wordlayout(x2$log2FoldChange,-log10(x2$padj),rownames(x2),cex=1)
text(nc[,1],nc[,2],label=rownames(x2),cex=.7)
legend("left", c("HR","NHEJ"),fill=c("red","blue"), bty="n") 
dev.off()

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


pdf("4_CD41-_untr_VS_CD41-_tr/volcano_labels.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41- Untreated / CD41- Treated )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[rownames(dds_res) %in% HR],
       -log10(dds_res$padj)[rownames(dds_res) %in% HR],
      col="red",pch=20)

points(dds_res$log2FoldChange[rownames(dds_res) %in% NHEJ],
       -log10(dds_res$padj)[rownames(dds_res) %in% NHEJ],
      col="blue",pch=20)
legend("topright", paste("CD41- Untreated:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41- Treated:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n") 

x=dds_res
x1=x[which((rownames(dds_res) %in% NHEJ) & abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05), ]
x2=x[which((rownames(dds_res) %in% HR) & abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05), ]
nc=wordlayout(x1$log2FoldChange,-log10(x1$padj),rownames(x1),cex=1)
text(nc[,1],nc[,2],label=rownames(x1),cex=.7)
nc=wordlayout(x2$log2FoldChange,-log10(x2$padj),rownames(x2),cex=1)
text(nc[,1],nc[,2],label=rownames(x2),cex=.7)
legend("left", c("HR","NHEJ"),fill=c("red","blue"), bty="n") 
dev.off()

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

pdf("5_CD41+_VS_CD41-_ControlledByTreatment/volcano_labels.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( CD41+ / CD41- )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[rownames(dds_res) %in% HR],
       -log10(dds_res$padj)[rownames(dds_res) %in% HR],
      col="red",pch=20)

points(dds_res$log2FoldChange[rownames(dds_res) %in% NHEJ],
       -log10(dds_res$padj)[rownames(dds_res) %in% NHEJ],
      col="blue",pch=20)
legend("topright", paste("CD41+:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("CD41-:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n") 

x=dds_res
x1=x[which((rownames(dds_res) %in% NHEJ) & abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05), ]
x2=x[which((rownames(dds_res) %in% HR) & abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05), ]
nc=wordlayout(x1$log2FoldChange,-log10(x1$padj),rownames(x1),cex=1)
text(nc[,1],nc[,2],label=rownames(x1),cex=.7)
nc=wordlayout(x2$log2FoldChange,-log10(x2$padj),rownames(x2),cex=1)
text(nc[,1],nc[,2],label=rownames(x2),cex=.7)
legend("left", c("HR","NHEJ"),fill=c("red","blue"), bty="n") 
dev.off()

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

pdf("6_untr_VS_tr_ControlledByCD41Status/volcano_labels.pdf")
plot(dds_res$log2FoldChange,-log10(dds_res$padj),xlab=expression('Log'[2]*' Fold Change ( Ctrl / Thpo )'),
              ylab=expression('-Log'[10]*' Q-values'),col=alpha("grey",.5),pch=20)
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.05),lty = 2,col="grey")
points(dds_res$log2FoldChange[rownames(dds_res) %in% HR],
       -log10(dds_res$padj)[rownames(dds_res) %in% HR],
      col="red",pch=20)

points(dds_res$log2FoldChange[rownames(dds_res) %in% NHEJ],
       -log10(dds_res$padj)[rownames(dds_res) %in% NHEJ],
      col="blue",pch=20)
legend("topright", paste("Ctrl:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 
legend("topleft", paste("Thpo:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.05))), bty="n") 

x=dds_res
x1=x[which((rownames(dds_res) %in% NHEJ) & abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05), ]
x2=x[which((rownames(dds_res) %in% HR) & abs(dds_res$log2FoldChange)>1 & dds_res$padj<0.05), ]
nc=wordlayout(x1$log2FoldChange,-log10(x1$padj),rownames(x1),cex=1)
text(nc[,1],nc[,2],label=rownames(x1),cex=.7)
nc=wordlayout(x2$log2FoldChange,-log10(x2$padj),rownames(x2),cex=1)
text(nc[,1],nc[,2],label=rownames(x2),cex=.7)
legend("left", c("HR","NHEJ"),fill=c("red","blue"), bty="n") 
dev.off()
