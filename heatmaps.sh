
cut -f 1,2,3 CD41+_over_CD41-.bed CD41+_tr_over_CD41+_untr.bed > CD412.bed
echo "# CD41+_tr specific" >> CD412.bed
cut -f 1,2,3 CD41+_untr_over_CD41+_tr.bed >> CD412.bed
echo "# CD41+_untr_over_CD41+_tr" >> CD412.bed
cut -f 1,2,3 CD41-_over_CD41+.bed >> CD412.bed
echo "# CD41-_over_CD41+ " >> CD412.bed
#cut -f 1,2,3 CD41_merge_NOT.bed >> CD41.bed
#echo "# common " >> CD41.bed


computeMatrix reference-point \
-S /root/ayako/ayako_dejavu/bw/CD41+_tr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41+_tr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_3.bw \
/root/ayako/ayako_dejavu/bw/CD41-_tr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41-_tr_2.bw \
-R CD412.bed --referencePoint center \
--sortRegions descend -bs 20 -a 1000 -b 1000 -p 40 -out CD412.mat

plotHeatmap --xAxisLabel "" --refPointLabel "ATAC Peak" --colorMap Blues -m CD412.mat -out CD412.pdf

##
##
