# Creating star index for RNA-Seq
~/myPrograms/STAR/bin/Linux_x86_64/STAR --runThreadN 50 --runMode genomeGenerate --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ \
--genomeFastaFiles /root/ayako/ref/mm10.fa \
--sjdbGTFfile /root/ayako/ref/gencode.vM15.annotation.gtf --sjdbOverhang 100
##
#
#################################
~/myPrograms/STAR/bin/Linux_x86_64/STAR --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ \
--readFilesCommand zcat \
--runThreadN 40 \
--readFilesIn /root/ayako/rna/Ctrl_CD41minus1_1.fastq.gz \
/root/ayako/rna/Ctrl_CD41minus1_2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ayako/bam/Ctrl_CD41minus1_

~/myPrograms/STAR/bin/Linux_x86_64/STAR --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ \
--readFilesCommand zcat \
--runThreadN 40 \
--readFilesIn /root/ayako/rna/Ctrl_CD41minus2_1.fastq.gz \
/root/ayako/rna/Ctrl_CD41minus2_2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ayako/bam/Ctrl_CD41minus2_

~/myPrograms/STAR/bin/Linux_x86_64/STAR --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ \
--readFilesCommand zcat \
--runThreadN 40 \
--readFilesIn /root/ayako/rna/Ctrl_CD41minus3_1.fastq.gz \
/root/ayako/rna/Ctrl_CD41minus3_2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ayako/bam/Ctrl_CD41minus3_
#############################################
~/myPrograms/STAR/bin/Linux_x86_64/STAR --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ \
--readFilesCommand zcat \
--runThreadN 40 \
--readFilesIn /root/ayako/rna/Ctrl_CD41plus1_1.fastq.gz \
/root/ayako/rna/Ctrl_CD41plus1_2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ayako/bam/Ctrl_CD41plus1_

~/myPrograms/STAR/bin/Linux_x86_64/STAR --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ \
--readFilesCommand zcat \
--runThreadN 40 \
--readFilesIn /root/ayako/rna/Ctrl_CD41plus2_1.fastq.gz \
/root/ayako/rna/Ctrl_CD41plus2_2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ayako/bam/Ctrl_CD41plus2_

~/myPrograms/STAR/bin/Linux_x86_64/STAR --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ \
--readFilesCommand zcat \
--runThreadN 40 \
--readFilesIn /root/ayako/rna/Ctrl_CD41plus3_1.fastq.gz \
/root/ayako/rna/Ctrl_CD41plus3_2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ayako/bam/Ctrl_CD41plus3_
############################################
~/myPrograms/STAR/bin/Linux_x86_64/STAR --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ \
--readFilesCommand zcat \
--runThreadN 40 \
--readFilesIn /root/ayako/rna/Ctrl_neg_4_1.fastq.gz \
/root/ayako/rna/Ctrl_neg_4_2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ayako/bam/Ctrl_neg_

###############################################

~/myPrograms/STAR/bin/Linux_x86_64/STAR --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ \
--readFilesCommand zcat \
--runThreadN 40 \
--readFilesIn /root/ayako/rna/Thpo_CD41minus1_1.fastq.gz \
/root/ayako/rna/Thpo_CD41minus1_2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ayako/bam/Thpo_CD41minus1_

~/myPrograms/STAR/bin/Linux_x86_64/STAR --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ \
--readFilesCommand zcat \
--runThreadN 40 \
--readFilesIn /root/ayako/rna/Thpo_CD41minus2_1.fastq.gz \
/root/ayako/rna/Thpo_CD41minus2_2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ayako/bam/Thpo_CD41minus2_

~/myPrograms/STAR/bin/Linux_x86_64/STAR --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ \
--readFilesCommand zcat \
--runThreadN 40 \
--readFilesIn /root/ayako/rna/Thpo_CD41minus3_1.fastq.gz \
/root/ayako/rna/Thpo_CD41minus3_2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ayako/bam/Thpo_CD41minus3_

###############################################

~/myPrograms/STAR/bin/Linux_x86_64/STAR --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ \
--readFilesCommand zcat \
--runThreadN 40 \
--readFilesIn /root/ayako/rna/Thpo_CD41plus1_1.fastq.gz \
/root/ayako/rna/Thpo_CD41plus1_2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ayako/bam/Thpo_CD41plus1_

~/myPrograms/STAR/bin/Linux_x86_64/STAR --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ \
--readFilesCommand zcat \
--runThreadN 40 \
--readFilesIn /root/ayako/rna/Thpo_CD41plus2_1.fastq.gz \
/root/ayako/rna/Thpo_CD41plus2_2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ayako/bam/Thpo_CD41plus2_

~/myPrograms/STAR/bin/Linux_x86_64/STAR --genomeDir /root/ayako/ref/mm10_star_sjdbO100/ \
--readFilesCommand zcat \
--runThreadN 40 \
--readFilesIn /root/ayako/rna/Thpo_CD41plus3_1.fastq.gz \
/root/ayako/rna/Thpo_CD41plus3_2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /root/ayako/bam/Thpo_CD41plus3_