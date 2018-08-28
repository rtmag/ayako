more gencode.vM15.annotation.gtf|grep -w "gene"|perl -pe 's/(.+)gene\_id.+gene\_name \"(.+)\"\; level.+/$1$2/g'| \
awk -F"\t" '{if($7=="+"){print $1"\t"$4"\t"$4+1"\t"$9"\t"$8"\t"$7};if($7=="-"){print $1"\t"$5-1"\t"$5"\t"$9"\t"$8"\t"$7}}' \
> gencode.vM15_tss.bed

############################################################################################################################
bed = read.table("gencode.vM15_tss.bed",sep="\t",stringsAsFactors=F)

x = read.table("../bulk_rna/gene_list_rna.txt",sep="\t",stringsAsFactors=F)
x = x[,1]

bed_m = bed[match(x, bed[,4]),]
write.table(bed_m,"gencode.vM15_tss_RNASEQ.bed",sep="\t",quote=F,col.names=F,row.names=F)
############################################################################################################################
more /root/ayako/ref/gencode.vM15_tss_RNASEQ.bed|awk -F"\t" '{print $1"\t"$2-500"\t"$2+500}'|grep -v "-" > m15_1kb_aroundTSS.bed
more /root/ayako/ref/gencode.vM15_tss_RNASEQ.bed|awk -F"\t" '{print $1"\t"$2-1000"\t"$2+1000}'|grep -v "-" > m15_2kb_aroundTSS.bed
more /root/ayako/ref/gencode.vM15_tss_RNASEQ.bed|awk -F"\t" '{print $1"\t"$2-2000"\t"$2+2000}'|grep -v "-" > m15_4kb_aroundTSS.bed


computeMatrix scale-regions \
-S /root/ayako/ayako_dejavu/bw/CD41+_tr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41+_tr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_3.bw \
/root/ayako/ayako_dejavu/bw/CD41-_tr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41-_tr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41-_untr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41-_untr_2.bw \
-R m15_1kb_aroundTSS.bed  \
--sortRegions descend -bs 1 -m 1 -p max -out m15_1kb_aroundTSS.mat --outFileNameMatrix m15_1kb_aroundTSS.txt

computeMatrix scale-regions \
-S /root/ayako/ayako_dejavu/bw/CD41+_tr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41+_tr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_3.bw \
/root/ayako/ayako_dejavu/bw/CD41-_tr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41-_tr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41-_untr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41-_untr_2.bw \
-R m15_2kb_aroundTSS.bed  \
--sortRegions descend -bs 1 -m 1 -p max -out m15_2kb_aroundTSS.mat --outFileNameMatrix m15_2kb_aroundTSS.txt

computeMatrix scale-regions \
-S /root/ayako/ayako_dejavu/bw/CD41+_tr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41+_tr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_3.bw \
/root/ayako/ayako_dejavu/bw/CD41-_tr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41-_tr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41-_untr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41-_untr_2.bw \
-R m15_4kb_aroundTSS.bed  \
--sortRegions descend -bs 1 -m 1 -p max -out m15_4kb_aroundTSS.mat --outFileNameMatrix m15_4kb_aroundTSS.txt
