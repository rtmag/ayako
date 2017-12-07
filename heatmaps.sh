
cut -f 1,2,3 CD41+_over_CD41-.bed CD41+_tr_over_CD41+_untr.bed > CD41.bed
echo "# CD41+_tr specific" >> CD41.bed
cut -f 1,2,3 CD41+_untr_over_CD41+_tr.bed >> CD41.bed
echo "# CD41+_untr_over_CD41+_tr" >> CD41.bed
cut -f 1,2,3 CD41-_over_CD41+.bed >> CD41.bed
echo "# CD41-_over_CD41+ " >> CD41.bed
cut -f 1,2,3 CD41_merge_NOT.bed >> CD41.bed
echo "# common " >> CD41.bed


computeMatrix reference-point \
-S /root/ayako/ayako_dejavu/bw/CD41+_tr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41+_tr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_3.bw \
/root/ayako/ayako_dejavu/bw/CD41-_tr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41-_tr_2.bw \
-R CD41.bed --referencePoint center \
--sortRegions descend -bs 20 -a 1000 -b 1000 -p 40 -out CD41.mat

plotHeatmap --xAxisLabel "" --refPointLabel "ATAC Peak" --colorMap Blues -m CD41.mat -out CD41.pdf

##
##
