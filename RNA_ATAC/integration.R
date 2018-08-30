library(graphics)

rna = readRDS("ayako_rpkm.rds")

rna_neg_untr = asinh(rowMeans(rna[,1:3]))
rna_plus_untr = asinh(rowMeans(rna[,4:6]))
rna_neg_tr = asinh(rowMeans(rna[,7:9]))
rna_plus_tr = asinh(rowMeans(rna[,10:12]))

pdf("ATAC_RNA_1kb_aroundTSS.pdf")
atac = read.table(pipe("more m15_1kb_aroundTSS.txt |grep -v '#'|perl -pe 's/genes:52459\t//g' "),sep="\t",header=T)
colnames(atac) = c("plus_tr_1","plus_tr_2","plus_untr_1","plus_untr_2","plus_untr_3","neg_tr_1","neg_tr_2","neg_untr_1","neg_untr_2")

par(mfrow=c(3,3))
# plus_tr_1
smoothScatter(asinh(atac[,1]),rna_plus_tr,nrpoints=0,xlab="RNA-Seq CD41+ treated",ylab="ATAC-Seq TSS CD41+ treated 1")
# plus_tr_2
smoothScatter(asinh(atac[,2]),rna_plus_tr,nrpoints=0,xlab="RNA-Seq CD41+ treated",ylab="ATAC-Seq TSS CD41+ treated 2")
# plus_untr_1
smoothScatter(asinh(atac[,3]),rna_plus_untr,nrpoints=0,xlab="RNA-Seq CD41+ untreated",ylab="ATAC-Seq TSS CD41+ untreated 1")
# plus_untr_2
smoothScatter(asinh(atac[,4]),rna_plus_untr,nrpoints=0,xlab="RNA-Seq CD41+ untreated",ylab="ATAC-Seq TSS CD41+ untreated 2")
# plus_untr_3
smoothScatter(asinh(atac[,5]),rna_plus_untr,nrpoints=0,xlab="RNA-Seq CD41+ untreated",ylab="ATAC-Seq TSS CD41+ untreated 3")
# neg_tr_1
smoothScatter(asinh(atac[,6]),rna_neg_tr,nrpoints=0,xlab="RNA-Seq CD41- treated",ylab="ATAC-Seq TSS CD41- treated 1")
# neg_tr_2
smoothScatter(asinh(atac[,7]),rna_neg_tr,nrpoints=0,xlab="RNA-Seq CD41- treated",ylab="ATAC-Seq TSS CD41- treated 2")
# neg_untr_1
smoothScatter(asinh(atac[,8]),rna_neg_untr,nrpoints=0,xlab="RNA-Seq CD41- untreated",ylab="ATAC-Seq TSS CD41- untreated 1")
# neg_untr_2
smoothScatter(asinh(atac[,9]),rna_neg_untr,nrpoints=0,xlab="RNA-Seq CD41- untreated",ylab="ATAC-Seq TSS CD41- untreated 2")
dev.off()

pdf("ATAC_RNA_2kb_aroundTSS.pdf")
atac = read.table(pipe("more m15_2kb_aroundTSS.txt |grep -v '#'|perl -pe 's/genes:52459\t//g' "),sep="\t",header=T)
colnames(atac) = c("plus_tr_1","plus_tr_2","plus_untr_1","plus_untr_2","plus_untr_3","neg_tr_1","neg_tr_2","neg_untr_1","neg_untr_2")

par(mfrow=c(3,3))
# plus_tr_1
smoothScatter(asinh(atac[,1]),rna_plus_tr,nrpoints=0,xlab="RNA-Seq CD41+ treated",ylab="ATAC-Seq TSS CD41+ treated 1")
# plus_tr_2
smoothScatter(asinh(atac[,2]),rna_plus_tr,nrpoints=0,xlab="RNA-Seq CD41+ treated",ylab="ATAC-Seq TSS CD41+ treated 2")
# plus_untr_1
smoothScatter(asinh(atac[,3]),rna_plus_untr,nrpoints=0,xlab="RNA-Seq CD41+ untreated",ylab="ATAC-Seq TSS CD41+ untreated 1")
# plus_untr_2
smoothScatter(asinh(atac[,4]),rna_plus_untr,nrpoints=0,xlab="RNA-Seq CD41+ untreated",ylab="ATAC-Seq TSS CD41+ untreated 2")
# plus_untr_3
smoothScatter(asinh(atac[,5]),rna_plus_untr,nrpoints=0,xlab="RNA-Seq CD41+ untreated",ylab="ATAC-Seq TSS CD41+ untreated 3")
# neg_tr_1
smoothScatter(asinh(atac[,6]),rna_neg_tr,nrpoints=0,xlab="RNA-Seq CD41- treated",ylab="ATAC-Seq TSS CD41- treated 1")
# neg_tr_2
smoothScatter(asinh(atac[,7]),rna_neg_tr,nrpoints=0,xlab="RNA-Seq CD41- treated",ylab="ATAC-Seq TSS CD41- treated 2")
# neg_untr_1
smoothScatter(asinh(atac[,8]),rna_neg_untr,nrpoints=0,xlab="RNA-Seq CD41- untreated",ylab="ATAC-Seq TSS CD41- untreated 1")
# neg_untr_2
smoothScatter(asinh(atac[,9]),rna_neg_untr,nrpoints=0,xlab="RNA-Seq CD41- untreated",ylab="ATAC-Seq TSS CD41- untreated 2")
dev.off()

pdf("ATAC_RNA_4kb_aroundTSS.pdf")
atac = read.table(pipe("more m15_4kb_aroundTSS.txt |grep -v '#'|perl -pe 's/genes:52459\t//g' "),sep="\t",header=T)
colnames(atac) = c("plus_tr_1","plus_tr_2","plus_untr_1","plus_untr_2","plus_untr_3","neg_tr_1","neg_tr_2","neg_untr_1","neg_untr_2")

par(mfrow=c(3,3))
# plus_tr_1
smoothScatter(asinh(atac[,1]),rna_plus_tr,nrpoints=0,xlab="RNA-Seq CD41+ treated",ylab="ATAC-Seq TSS CD41+ treated 1")
# plus_tr_2
smoothScatter(asinh(atac[,2]),rna_plus_tr,nrpoints=0,xlab="RNA-Seq CD41+ treated",ylab="ATAC-Seq TSS CD41+ treated 2")
# plus_untr_1
smoothScatter(asinh(atac[,3]),rna_plus_untr,nrpoints=0,xlab="RNA-Seq CD41+ untreated",ylab="ATAC-Seq TSS CD41+ untreated 1")
# plus_untr_2
smoothScatter(asinh(atac[,4]),rna_plus_untr,nrpoints=0,xlab="RNA-Seq CD41+ untreated",ylab="ATAC-Seq TSS CD41+ untreated 2")
# plus_untr_3
smoothScatter(asinh(atac[,5]),rna_plus_untr,nrpoints=0,xlab="RNA-Seq CD41+ untreated",ylab="ATAC-Seq TSS CD41+ untreated 3")
# neg_tr_1
smoothScatter(asinh(atac[,6]),rna_neg_tr,nrpoints=0,xlab="RNA-Seq CD41- treated",ylab="ATAC-Seq TSS CD41- treated 1")
# neg_tr_2
smoothScatter(asinh(atac[,7]),rna_neg_tr,nrpoints=0,xlab="RNA-Seq CD41- treated",ylab="ATAC-Seq TSS CD41- treated 2")
# neg_untr_1
smoothScatter(asinh(atac[,8]),rna_neg_untr,nrpoints=0,xlab="RNA-Seq CD41- untreated",ylab="ATAC-Seq TSS CD41- untreated 1")
# neg_untr_2
smoothScatter(asinh(atac[,9]),rna_neg_untr,nrpoints=0,xlab="RNA-Seq CD41- untreated",ylab="ATAC-Seq TSS CD41- untreated 2")
dev.off()

#

pdf("ATAC_RNA_1kb_aroundTSS.pdf")
atac = read.table(pipe("more m15_1kb_aroundTSS_sum.txt |grep -v '#'|perl -pe 's/genes:52459\t//g' "),sep="\t",header=T)
colnames(atac) = c("plus_tr_1","plus_tr_2","plus_untr_1","plus_untr_2","plus_untr_3","neg_tr_1","neg_tr_2","neg_untr_1","neg_untr_2")

par(mfrow=c(3,3))
# plus_tr_1
smoothScatter(rna_plus_tr,asinh(atac[,1]),nrpoints=0,xlab="RNA-Seq CD41+ treated",ylab="ATAC-Seq TSS CD41+ treated 1")
# plus_tr_2
smoothScatter(rna_plus_tr,asinh(atac[,2]),nrpoints=0,xlab="RNA-Seq CD41+ treated",ylab="ATAC-Seq TSS CD41+ treated 2")
# plus_untr_1
smoothScatter(rna_plus_untr,asinh(atac[,3]),nrpoints=0,xlab="RNA-Seq CD41+ untreated",ylab="ATAC-Seq TSS CD41+ untreated 1")
# plus_untr_2
smoothScatter(rna_plus_untr,asinh(atac[,4]),nrpoints=0,xlab="RNA-Seq CD41+ untreated",ylab="ATAC-Seq TSS CD41+ untreated 2")
# plus_untr_3
smoothScatter(rna_plus_untr,asinh(atac[,5]),nrpoints=0,xlab="RNA-Seq CD41+ untreated",ylab="ATAC-Seq TSS CD41+ untreated 3")
# neg_tr_1
smoothScatter(rna_neg_tr,asinh(atac[,6]),nrpoints=0,xlab="RNA-Seq CD41- treated",ylab="ATAC-Seq TSS CD41- treated 1")
# neg_tr_2
smoothScatter(rna_neg_tr,asinh(atac[,7]),nrpoints=0,xlab="RNA-Seq CD41- treated",ylab="ATAC-Seq TSS CD41- treated 2")
# neg_untr_1
smoothScatter(rna_neg_untr,asinh(atac[,8]),nrpoints=0,xlab="RNA-Seq CD41- untreated",ylab="ATAC-Seq TSS CD41- untreated 1")
# neg_untr_2
smoothScatter(rna_neg_untr,asinh(atac[,9]),nrpoints=0,xlab="RNA-Seq CD41- untreated",ylab="ATAC-Seq TSS CD41- untreated 2")
dev.off()


smoothScatter(rna_neg_untr[rna_neg_untr>1],asinh(atac[rna_neg_untr>1,9]),nrpoints=0,xlab="RNA-Seq CD41- untreated",ylab="ATAC-Seq TSS CD41- untreated 2")

