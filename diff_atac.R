cat CD41+_tr_fdr5_peaks.broadPeak CD41+_untr_fdr5_peaks.broadPeak CD41-_tr_fdr5_peaks.broadPeak| \
sort -k1,1 -k2,2n |bedtools merge -i - > CD41_merged_peaks.bed

#

bed_to_granges <- function(file){
   df <- read.table(file,
                    header=F,
                    stringsAsFactors=F)
 
   if(length(df) > 6){
      df <- df[,-c(7:length(df))]
   }
 
   if(length(df)<3){
      stop("File has less than 3 columns")
   }
 
   header <- c('chr','start','end','id','score','strand')
   names(df) <- header[1:length(names(df))]
 
   if('strand' %in% colnames(df)){
      df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
   }
 
   library("GenomicRanges")
 
   if(length(df)==3){
      gr <- with(df, GRanges(chr, IRanges(start, end)))
   } else if (length(df)==4){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
   } else if (length(df)==5){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
   } else if (length(df)==6){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
   }
   return(gr)
}

#regions=bed_to_granges("BROADPEAK_TSS.bed")
regions=bed_to_granges("BROADPEAK_66.bed")
 counts <- regionCounts(bam.files, regions, ext=300, param=param)
countData=assay(counts)
colnames(countData)=bam.files

require(DESeq2)

colData <- data.frame(group=colnames(countData))
dds <- DESeqDataSetFromMatrix(
       countData = countData,
       colData = colData,
       design = ~ group)
dds <- estimateSizeFactors(dds)
sizeFactors(dds) <- normfacs


rld <- rlogTransformation( dds )

x=assay(rld)
postscript("~/Downloads/h2az/66_acH2AZ_H2AZ.ps")
boxplot((x[,3]/x[,1]),(x[,4]/x[,1]),outline=0,border=c("blue","goldenrod4"),width=c(4,4),names=c("siControl","siTIP60A"),ylab="acH2Az/H2AZ across peaks at TSS",boxlwd = 3)
dev.off()

