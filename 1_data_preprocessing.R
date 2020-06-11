#Thie is the preprocessing step.
#TAD boundaries were identified from Hi-C data
#Narrow peaks were extracted from ChIPseq data
#Transcription factor binding sites were extracted from wgEncodeRegTfbsClusteredV3.bed file
#Transcription start site information for coding, noncoding, and housekeeping genes were extracted from tss_RefGene.txt





library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(data.table)
library(stringr)
library(plyr)

#Get TAD boundaries from TAD file
TAD_T47D<- read_table2("Data/ENCODE3-T470-HindIII__hg19__genome__C-40000-iced.tads.bed",col_names = F)
bdsdf_T47D <-matrix("0",nrow=3189,ncol=3)
colnames(bdsdf_T47D)<-c("seqnames","start","end")
for(i in 1:(nrow(TAD_T47D)-1)){
  if(TAD_T47D[i,"X1"] == TAD_T47D[i+1,"X1"])
  {
    bdsdf_T47D[i,"seqnames"]<-as.character(TAD_T47D[i,"X1"])
    bdsdf_T47D[i,"start"]<-(as.numeric(TAD_T47D[i,"X3"])+as.numeric(TAD_T47D[i+1,"X2"])+1)/2 - 100000 #boundary width is 200k
    bdsdf_T47D[i,"end"]<-(as.numeric(TAD_T47D[i,"X3"])+as.numeric(TAD_T47D[i+1,"X2"])+1)/2 + 100000
  }
}

bdsdf_T47D<-bdsdf_T47D[-(which(bdsdf_T47D[,"seqnames"]=="0")),] #remove the rows inbetween chromosomes

bdsGR_T47D <- GRanges(seqnames = bdsdf_T47D[,"seqnames"],
                      ranges = IRanges(start=as.numeric(bdsdf_T47D[,"start"]),end=as.numeric(bdsdf_T47D[,"end"])),
                      strand='*')


####################################################################################
#########         beds files      
####################################################################################


bedfl<-read_table2('Data/T47D_hg19/List-of-narrowPeaks-T47D.txt',col_names = T)
bed_data<-lapply(1:dim(bedfl)[1]
                 ,function(i){read_table2(paste('Data/T47D_hg19/',bedfl$file_name[i],'.bed.gz',sep=''),col_names = F)}
)


bedf2<-read_table2('Data/T47D_hg19/List-of-bam-T47D.txt',col_names = T)
bed_data2<-lapply(1:dim(bedf2)[1]
                  ,function(i){read_table2(paste('Data/T47D_hg19/bam_files/macs2/',bedf2$file_name[i],'_peaks.narrowPeak',sep=''),col_names = F)}
)


bedf3<-read_table2('Data/T47D_hg38/List-of-narrowPeaks-hg38.txt',col_names = T)
chain<-import.chain('Data/Annotations/liftOver/hg38ToHg19.over.chain')
for (i in 1:dim(bedf3)[1]){
  x1 <- read.table(paste("Data/T47D_hg38/",bedf3$file_name[i],".bed.gz",sep=''),header=F)
  x2 <- GRanges(seqnames=x1$V1,
                ranges=IRanges(start=x1$V2,end=x1$V3),
                strand='*',
                V4=x1$V4,
                V5=x1$V5,
                V6=x1$V6,
                V7=x1$V7,
                V8=x1$V8)
  x3 <- liftOver(x2,chain)
  x4<-data.frame(unlist(x3),stringsAsFactors = F)
  write.table(x4,file=paste("Data/ChIPseq-out/T47D_hg38/liftover-hg38-to-hg19/",bedf3$file_name[i],".hg19.bed",sep=''),quote=F,sep='\t',col.names=F, row.names=F)
}
rm(x1,x2,x3,x4)
bed_data3<-lapply(1:dim(bedf3)[1]
                  ,function(i){read_table2(paste('Data/ChIPseq-out/T47D_hg38/liftover-hg38-to-hg19/',bedf3$file_name[i],'.hg19.bed',sep=''),col_names = F)})

peak_col_names<-c('chrom','start','end','name','score','strand','sig','pv','qv','peak')

###############################################
###     TFBS information (Transcription factors binding sites)
##############################################


BS<-fread('R_codes_ZZG_for_TADs/0821/mcf7_chip/BindingSites/wgEncodeRegTfbsClusteredV3.bed',sep = '\t',header = F)
BSwC<-fread('../R_codes_ZZG_for_TADs/0821/mcf7_chip/BindingSites/wgEncodeRegTfbsClusteredWithCellsV3.bed',sep = '\t',header = F)
names(BSwC)<-c('chr','start','end','type','value','cell')
names(BS)<-c('chr','start','end','type','value','NoC','cellNo1','cellNo2')

tfbsGR<-GRanges(seqnames=BSwC$chr
                ,ranges=IRanges(start=BSwC$start,end=BSwC$end)
                ,strand='*'
                ,type=BSwC$type
                ,value=BSwC$value)

tfbsname<-unique(tfbsGR$type)


####################################################################
#####          TSS information
####################################################################
tss<-read_tsv('../R_codes_ZZG_for_TADs/0821/tss_RefGene.txt')
hks<-read_table2('../R_codes_ZZG_for_TADs/0821/tss_HK_genes.txt',col_names = c('name','id'))
tss<-tss[tss$chrom %in% chrnames,]

#using the  unique row in tssfile 
tmp<-duplicated(tss[,c('chrom','txStart','txEnd','name2')])
tss<-tss[!tmp,]

#location of non-coding Tss
loc_non_coding<-str_detect(tss$name2,'LINC*|LOC*|MIR*')

#location of coding Tss
loc_coding<-!loc_non_coding

#location of HK Tss
loc_hk<-tss$name2 %in% unique(hks$name)

####deal with non_coding Tss with bins
tssGR<-GRanges(seqnames=tss$chrom,
               ranges=IRanges(start=tss$txStart,end=tss$txEnd),
               strand=tss$strand)


###################################################
###   location info 
##################################################
##get center and chr_size information for each chrs
chr_info_1<-read_tsv('../R_codes_ZZG_for_TADs/0821/mcf7_chip/center_chr/center.txt')
chr_info_2<-read_tsv('../R_codes_ZZG_for_TADs/0821/mcf7_chip/center_chr/chr_size.txt',col_names = F)



chrnames<-c(paste('chr',1:22,sep=''),'chrX')

save(list=ls(),file='pre_info_12102019.Rdata')
