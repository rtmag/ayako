

# for file in /root/ayako/single_cell/giga_38e_cat/*; do
#    trim_galore --quality 20 --illumina -o /root/ayako/single_cell/giga_38e_cat_trim/ $file
# done

# for file in /root/ayako/single_cell/giga_1ad_cat/*; do
#    trim_galore --quality 20 --illumina -o /root/ayako/single_cell/giga_1ad_cat_trim/ $file
# done

#############################################################

~/myPrograms/STAR/bin/Linux_x86_64/STAR --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ --genomeLoad LoadAndExit


for file in /root/ayako/single_cell/giga_38e_cat_trim/*_trimmed.fq.gz; do

    file_name=${file/_trimmed.fq.gz//}
    file_name=${file_name/\/root\/ayako\/single_cell\/giga_38e_cat_trim\///}
    file_name=${file_name////}
            
~/myPrograms/STAR/bin/Linux_x86_64/STAR --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ \
--genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 \
--readFilesCommand zcat \
--runThreadN 40 \
--readFilesIn $file \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5 \
--outFileNamePrefix /root/ayako/single_cell/giga_38e_bam\/$file_name\_
    
done


for file in /root/ayako/single_cell/giga_1ad_cat_trim/*_trimmed.fq.gz; do

    file_name=${file/_trimmed.fq.gz//}
    file_name=${file_name/\/root\/ayako\/single_cell\/giga_38e_cat_trim\///}
    file_name=${file_name////}
            
~/myPrograms/STAR/bin/Linux_x86_64/STAR --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ \
--genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 \
--readFilesCommand zcat \
--runThreadN 40 \
--readFilesIn $file \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5 \
--outFileNamePrefix /root/ayako/single_cell/giga_1ad_bam\/$file_name\_
    
done


~/myPrograms/STAR/bin/Linux_x86_64/STAR --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ --genomeLoad Remove
##
#




