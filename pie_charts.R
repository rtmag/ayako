pdf("pie_1_CD41+_untr_over_CD41+_tr-mm10_annotation.pdf")
res=read.table(pipe("more 1_CD41+_untr_over_CD41+_tr-mm10_annotation.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>2]
pie(sort(tdown), main=,cex=.8)
dev.off()

pdf("pie_1_CD41+_tr_over_CD41+_untr-mm10_annotation.pdf")
res=read.table(pipe("more 1_CD41+_tr_over_CD41+_untr-mm10_annotation.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>50]
pie(sort(tdown), main=,cex=.8)
dev.off()

pdf("pie_2_CD41-_tr_over_CD41-_untr-mm10_annotation.pdf")
res=read.table(pipe("more 2_CD41-_tr_over_CD41-_untr-mm10_annotation.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>20]
pie(sort(tdown), main=,cex=.8)
dev.off()

pdf("pie_2_CD41-_untr_over_CD41-_tr-mm10_annotation.pdf")
res=read.table(pipe("more 2_CD41-_untr_over_CD41-_tr-mm10_annotation.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>20]
pie(sort(tdown), main=,cex=.8)
dev.off()





pdf("pie_3_CD41-_untr_over_CD41+_untr-mm10_annotation.pdf")
res=read.table(pipe("more 3_CD41-_untr_over_CD41+_untr-mm10_annotation.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>2]
pie(sort(tdown), main=,cex=.8)
dev.off()

pdf("pie_3_CD41+_untr_over_CD41-_untr-mm10_annotation.pdf")
res=read.table(pipe("more 3_CD41+_untr_over_CD41-_untr-mm10_annotation.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>2]
pie(sort(tdown), main=,cex=.8)
dev.off()

pdf("pie_4_CD41-_tr_over_CD41+_tr-mm10_annotation.pdf")
res=read.table(pipe("more 4_CD41-_tr_over_CD41+_tr-mm10_annotation.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>50]
pie(sort(tdown), main=,cex=.8)
dev.off()

pdf("pie_4_CD41+_tr_over_CD41-_tr-mm10_annotation.pdf")
res=read.table(pipe("more 4_CD41+_tr_over_CD41-_tr-mm10_annotation.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>2]
pie(sort(tdown), main=,cex=.8)
dev.off()
