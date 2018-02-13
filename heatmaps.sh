
#cut -f 1,2,3 CD41+_over_CD41-.bed CD41+_tr_over_CD41+_untr.bed > CD412.bed
#echo "# CD41+_tr specific" >> CD412.bed
#cut -f 1,2,3 CD41+_untr_over_CD41+_tr.bed >> CD412.bed
#echo "# CD41+_untr_over_CD41+_tr" >> CD412.bed
#cut -f 1,2,3 CD41-_over_CD41+.bed >> CD412.bed
#echo "# CD41-_over_CD41+ " >> CD412.bed
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
/root/ayako/ayako_dejavu/bw/CD41-_untr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41-_untr_2.bw \
-R /root/ayako/ayako_dejavu/diffbind/diffbind_deseq2/ATAC-Seq_merged_LRT_FDR5_.bed --referencePoint center \
--sortRegions descend -bs 20 -a 1000 -b 1000 -p 40 -out /root/ayako/ayako_dejavu/heatmaps/ATAC-Seq_merged_LRT_FDR5.mat

plotHeatmap --xAxisLabel "" --refPointLabel "ATAC Peak" --colorMap Blues -m /root/ayako/ayako_dejavu/heatmaps/ATAC-Seq_merged_LRT_FDR5.mat -out /root/ayako/ayako_dejavu/heatmaps/ATAC-Seq_merged_LRT_FDR5.pdf



computeMatrix reference-point \
-S /root/ayako/ayako_dejavu/bw/CD41+_tr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41+_tr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_3.bw \
/root/ayako/ayako_dejavu/bw/CD41-_tr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41-_tr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41-_untr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41-_untr_2.bw \
-R /root/ayako/ayako_dejavu/diffbind/diffbind_deseq2/ATAC-Seq_merged_LRT_FDR1_.bed --referencePoint center \
--sortRegions descend -bs 20 -a 1000 -b 1000 -p 40 -out /root/ayako/ayako_dejavu/heatmaps/ATAC-Seq_merged_LRT_FDR1.mat

plotHeatmap --xAxisLabel "" --refPointLabel "ATAC Peak" --colorMap Blues -m /root/ayako/ayako_dejavu/heatmaps/ATAC-Seq_merged_LRT_FDR1.mat -out /root/ayako/ayako_dejavu/heatmaps/ATAC-Seq_merged_LRT_FDR1.pdf

##
##
