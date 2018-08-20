
countData=readRDS('atac_countdata.rds')

colnames(countData)=c("CD41_plus_untr_1","CD41_plus_untr_2","CD41_plus_untr_3",
                      "CD41_plus_tr_1","CD41_plus_tr_2",
                      "CD41_minus_tr_1","CD41_minus_tr_2",
                      "CD41_minus_untr_1","CD41_minus_untr_2")
require(DESeq2)


design<-data.frame(CD41=c("CD41plus","CD41plus","CD41plus",
                           "CD41plus","CD41plus",
                           "CD41minus","CD41minus",
                           "CD41minus","CD41minus"),
                   Treatment=c("Ctrl","Ctrl","Ctrl",
                           "Thpo","Thpo",
                           "Thpo","Thpo",
                           "Ctrl","Ctrl") )

dds <- DESeqDataSetFromMatrix(countData = countData, colData = design, 
                  design = ~ CD41 + Treatment )

dds <- DESeq(dds, test="LRT", 
           full= ~ CD41 + Treatment, 
           reduced= ~ Treatment )
res <- results(dds,contrast=c("CD41","CD41plus","CD41minus"))

ix=res$padj<0.05 & res$log2FoldChange<(-1)
ix[is.na(ix)]=FALSE
write.table(gsub("_","\t",rownames(res[ix,])),"CD41-_over_CD41+.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)

ix=res$padj<0.05 & res$log2FoldChange>(1)
ix[is.na(ix)]=FALSE
write.table(gsub("_","\t",rownames(res[ix,])),"CD41+_over_CD41-.bed",quote=FALSE,col.names=FALSE,row.names=FALSE)

cut -f 1,3,5 CD41-_over_CD41+.bed > 5_CD41-_over_CD41+.bed 
cut -f 1,3,5 CD41+_over_CD41-.bed > 5_CD41+_over_CD41-.bed 
