# Creating star index for RNA-Seq
~/myPrograms/STAR/bin/Linux_x86_64/STAR --runThreadN 50 --runMode genomeGenerate --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ \
--genomeFastaFiles /root/ayako/ref/mm10.fa \
--sjdbGTFfile /root/ayako/ref/gencode.vM15.annotation.gtf --sjdbOverhang 100
##
#
