
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

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC Peak" --colorMap RdBu_r --kmeans 4 -m /root/ayako/ayako_dejavu/heatmaps/ATAC-Seq_merged_LRT_FDR5.mat -out /root/ayako/ayako_dejavu/heatmaps/ATAC-Seq_merged_LRT_FDR5.pdf



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

cat 1_CD41+_tr_over_CD41+_untr.bed > combined_CD41+_tr_vs_CD41+_untr.bed
echo "#CD41+ treated Open regions" >> combined_CD41+_tr_vs_CD41+_untr.bed
cat 1_CD41+_untr_over_CD41+_tr.bed >> combined_CD41+_tr_vs_CD41+_untr.bed
echo "#CD41+ untreated Open regions" >> combined_CD41+_tr_vs_CD41+_untr.bed

cat 2_CD41-_tr_over_CD41-_untr.bed > combined_CD41-_tr_vs_CD41-_untr.bed
echo "#CD41- treated Open regions" >> combined_CD41-_tr_vs_CD41-_untr.bed
cat 2_CD41-_untr_over_CD41-_tr.bed >> combined_CD41-_tr_vs_CD41-_untr.bed
echo "#CD41- untreated Open regions" >> combined_CD41-_tr_vs_CD41-_untr.bed

cat 3_CD41-_untr_over_CD41+_untr.bed > combined_CD41-_untr_vs_CD41+_untr.bed
echo "#CD41- untreated Open regions" >> combined_CD41-_untr_vs_CD41+_untr.bed
cat 3_CD41+_untr_over_CD41-_untr.bed >> combined_CD41-_untr_vs_CD41+_untr.bed
echo "#CD41+ untreated Open regions" >> combined_CD41-_untr_vs_CD41+_untr.bed

cat 4_CD41-_tr_over_CD41+_tr.bed > combined_CD41-_tr_vs_CD41+_tr.bed
echo "#CD41- treated Open regions" >> combined_CD41-_tr_vs_CD41+_tr.bed
cat 4_CD41+_tr_over_CD41-_tr.bed >> combined_CD41-_tr_vs_CD41+_tr.bed
echo "#CD41+ treated Open regions" >> combined_CD41-_tr_vs_CD41+_tr.bed

#####################################


computeMatrix reference-point \
-S /root/ayako/ayako_dejavu/bw/CD41+_tr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41+_tr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_3.bw \
-R combined_CD41+_tr_vs_CD41+_untr.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p 40 -out combined_CD41+_tr_vs_CD41+_untr.mat


computeMatrix reference-point \
-S \
/root/ayako/ayako_dejavu/bw/CD41-_tr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41-_tr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41-_untr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41-_untr_2.bw \
-R combined_CD41-_tr_vs_CD41-_untr.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p 40 -out combined_CD41-_tr_vs_CD41-_untr.mat


computeMatrix reference-point \
-S \
/root/ayako/ayako_dejavu/bw/CD41-_untr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41-_untr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41+_untr_3.bw \
-R combined_CD41-_untr_vs_CD41+_untr.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p 40 -out combined_CD41-_untr_vs_CD41+_untr.mat


computeMatrix reference-point \
-S \
/root/ayako/ayako_dejavu/bw/CD41-_tr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41-_tr_2.bw \
/root/ayako/ayako_dejavu/bw/CD41+_tr_1.bw \
/root/ayako/ayako_dejavu/bw/CD41+_tr_2.bw \
-R combined_CD41-_tr_vs_CD41+_tr.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p 40 -out combined_CD41-_tr_vs_CD41+_tr.mat


plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC-Seq Peak" --colorMap RdBu_r \
-m combined_CD41+_tr_vs_CD41+_untr.mat \
 --samplesLabel "CD41+ treated 1" "CD41+ treated 2" "CD41+ untreated 1" "CD41+ untreated 2" "CD41+ untreated 3"  \
-out combined_CD41+_tr_vs_CD41+_untr.pdf

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC-Seq Peak" --colorMap RdBu_r \
-m combined_CD41-_tr_vs_CD41-_untr.mat \
 --samplesLabel "CD41- treated 1" "CD41- treated 2" "CD41- untreated 1" "CD41- untreated 2"  \
-out combined_CD41-_tr_vs_CD41-_untr.pdf

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC-Seq Peak" --colorMap RdBu_r \
-m combined_CD41-_untr_vs_CD41+_untr.mat \
 --samplesLabel "CD41- untreated 1" "CD41- untreated 2" "CD41+ untreated 1" "CD41+ untreated 2" "CD41+ untreated 3"  \
-out combined_CD41-_untr_vs_CD41+_untr.pdf

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "ATAC-Seq Peak" --colorMap RdBu_r \
-m combined_CD41-_tr_vs_CD41+_tr.mat \
 --samplesLabel "CD41- treated 1" "CD41- treated 2" "CD41+ treated 1" "CD41+ treated 2" \
-out combined_CD41-_tr_vs_CD41+_tr.pdf
#############
##
#
