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

# Boxplotter

boxploter=function(name_gene){
bdata= vsd[rownames(countData) %in% name_gene]
                                 
Ctrl_CD41minus = 1:3
Ctrl_CD41plus = 4:6
Thpo_CD41minus = 7:9
Thpo_CD41plus = 10:12
Ctrl_neg = 13
  
boxplot(bdata[Ctrl_CD41plus],bdata[Thpo_CD41plus],bdata[Ctrl_CD41minus],bdata[Thpo_CD41minus],bdata[Ctrl_neg],
        names=c("Ctrl_CD41plus","Thpo_CD41plus","Ctrl_CD41minus","Thpo_CD41minus","Ctrl_neg"), 
        main = name_gene, col=c("#ffb3ba","#baffc9","#bae1ff"),outline=F)
} 


pdf("boxplot_major_HSC_genes.pdf")
par(mfrow=c(3,3))
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

pdf("boxplot_cell_cycle_genes.pdf")
par(mfrow=c(3,3))
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

pdf("boxplot_major_MK_genes.pdf")
par(mfrow=c(3,3))
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
par(mfrow=c(3,3))
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

