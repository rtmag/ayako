

for file in /root/ayako/single_cell/giga_38e_cat/*; do
    trim_galore --quality 20 --illumina -o /root/ayako/single_cell/giga_38e_cat_trim/ $file
done

for file in /root/ayako/single_cell/giga_1ad_cat/*; do
    trim_galore --quality 20 --illumina -o /root/ayako/single_cell/giga_1ad_cat_trim/ $file
done
