library(graphics)

atac = read.table(pipe("more m15_1kb_aroundTSS.txt |grep -v '#'|perl -pe 's/genes:52457\t//g' "),sep="\t",header=T)
colnames(atac) = c("plus_tr_1","plus_tr_2","plus_untr_1","plus_untr_2","plus_untr_3","neg_tr_1","neg_tr_2","neg_untr_1","neg_untr_2")

rna = readRDS("ayako_rpkm.rds")

rna_neg_untr = asinh(rowMeans(rna[,1:3]))
rna_plus_untr = asinh(rowMeans(rna[,4:6]))
rna_neg_tr = asinh(rowMeans(rna[,7:9]))
rna_plus_tr = asinh(rowMeans(rna[,10:12]))

# [1] "plus_tr_1"   "plus_tr_2"   "plus_untr_1" "plus_untr_2" "plus_untr_3"
[6] "neg_tr_1"    "neg_tr_2"    "neg_untr_1"  "neg_untr_2" 

# plus_tr_1
smoothScatter(rna_plus_tr,atac[,1],nrpoints=0)

# plus_tr_2
rna_plus_tr,atac[,2]

# plus_untr_1
rna_plus_untr,atac[,3]

# plus_untr_2
rna_plus_untr,atac[,4]

# plus_untr_3
rna_plus_untr,atac[,5]

# neg_tr_1
rna_neg_tr,atac[,6]

# neg_tr_2
rna_neg_tr,atac[,7]

# neg_untr_1
rna_neg_untr,atac[,8]

# neg_untr_2
rna_neg_untr,atac[,9]
