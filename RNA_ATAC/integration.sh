more gencode.vM15.annotation.gtf|grep -w "gene"|perl -pe 's/(.+)gene\_id.+gene\_name \"(.+)\"\; level.+/$1$2/g'| \
awk -F"\t" '{if($7=="+"){print $1"\t"$4"\t"$4+1"\t"$9"\t"$8"\t"$7};if($7=="-"){print $1"\t"$5-1"\t"$5"\t"$9"\t"$8"\t"$7}}' \
> gencode.vM15_tss.bed

############################################################################################################################
bed = read.table("gencode.vM15_tss.bed",sep="\t",stringsAsFactors=F)

x = read.table("../bulk_rna/gene_list_rna.txt",sep="\t",stringsAsFactors=F)
x = x[,1]

bed_m = bed[match(x, bed[,4]),]
write.table(bed_m,"gencode.vM15_tss_RNASEQ.bed",sep="\t",quote=F,col.names=F,row.names=F)

