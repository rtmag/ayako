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

pdf("Diagnostic_pca.pdf")
plotPCA(dLRT_vsd,ntop=50000,intgroup=c("group"))
dev.off()


