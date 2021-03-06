more metadata.tsv |head -n2|cut -f1,7,8,10,12,13,38

####
mkdir lifton

# uncompressing and cut 1-3 bed
more metadata.tsv|cut -f1,7,8,10,12,13,38|grep "mm9" |cut -f1| perl -pe 's/(.+)\n/zcat $1.bed.gz \| cut \-f 1\-3 \> mm9_bed\/$1.bed\n/g'

# liftovering
ls -1 mm9_bed/| perl -pe 's/(.+)\n/\~\/myPrograms\/liftOver mm9\_bed\/$1 hg19ToHg38\.over\.chain mm9\_mm10\/$1 mm9\_mm10\/unknown\.bed \n/g'

#mm 10
more metadata.tsv|cut -f1,7,8,10,12,13,38|grep "mm10" |cut -f1| perl -pe 's/(.+)\n/zcat $1.bed.gz \| cut \-f 1\-3 \> mm9\_mm10\/$1.bed\n/g'

# anno table
more metadata.tsv|cut -f1,7,8,10,12,13|tail -n +2|head -n22

more metadata.tsv|cut -f1|tail -n +2 > file_table.txt
more metadata.tsv|cut -f7,8,10,12,13|tail -n +2|perl -pe 's/\t/\_/g'|perl -pe 's/\_+/\_/g' > anno_table.txt
paste file_table.txt anno_table.txt >anno.txt

while read line; do
   file=$(echo $line | cut -d " " -f 1)
   caca=$(echo $line | cut -d " " -f 2-999)
   echo "perl -pe 's/\n/\t$caca\n/g' mm9_mm10/$file.bed > mm10ano/$file.bed"
done < anno.txt > finalstep.bash

cat * > encode_chip-seq_mm10_unsorted.bed
sort -k1,1 -k2,2n encode_chip-seq_mm10_unsorted.bed > encode_chip-seq_mm10.bed

~/Downloads/bedToBigBed encode_chip-seq_mm10.bed mm10.chrom.sizes encode_chip-seq_mm10.bb
