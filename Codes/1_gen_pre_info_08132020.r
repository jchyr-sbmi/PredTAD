#GM12878 cell line
#with methylation RAD21, ATACseq 08132020


library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(data.table)
library(stringr)
library(plyr)



#Get TAD boundaries from TAD file
TAD<- read_table2("Data/Hi-C/GM12878_Lieberman-raw_TADs.txt",col_names = F)
#colnames(TAD)[1:3]<- c("X1","X2","X3") #chr, start, stop
bdsdf <-matrix("0",nrow=nrow(TAD),ncol=3)
colnames(bdsdf)<-c("seqnames","start","end")
for(i in 1:(nrow(TAD)-1)){
  if(TAD[i,"X1"] == TAD[i+1,"X1"])
    {
    bdsdf[i,"seqnames"]<-as.character(TAD[i,"X1"])
    bdsdf[i,"start"]<-(as.numeric(TAD[i,"X3"])+as.numeric(TAD[i+1,"X2"])+1)/2 - 100000 #boundary width is 200k
    bdsdf[i,"end"]<-(as.numeric(TAD[i,"X3"])+as.numeric(TAD[i+1,"X2"])+1)/2 + 100000
    }
}

bdsdf<-bdsdf[-(which(bdsdf[,"seqnames"]=="0")),] #remove the rows inbetween chromosomes

bdsGR <- GRanges(seqnames = bdsdf[,"seqnames"],
                      ranges = IRanges(start=as.numeric(bdsdf[,"start"]),end=as.numeric(bdsdf[,"end"])),
                      strand='*')
					  
bdsGR_GM12878 <- bdsGR

####################################################################################
#########         beds files      
####################################################################################

bedfl<-read_table2('Data/ENCODE/GM12878_List-of-Bed_narrowPeaks_hg19.txt',col_names = T)
bed_data1<-lapply(1:dim(bedfl)[1]
                 ,function(i){read_table2(paste('Data/ENCODE/',bedfl$file_name[i],'.bed.gz',sep=''),col_names = F)})
#standardize


peak_col_names<-c('chrom','start','end','name','score','strand','sig','pv','qv','peak')

###############################################
###     TFBS information (Transcription factors binding sites)
##############################################

BS<-fread('Data/TFBS/wgEncodeRegTfbsClusteredV3.bed',sep = '\t',header = F)
BSwC<-fread('Data/TFBS/wgEncodeRegTfbsClusteredWithCellsV3.bed',sep = '\t',header = F)
names(BSwC)<-c('chr','start','end','type','value','cell')
names(BS)<-c('chr','start','end','type','value','NoC','cellNo1','cellNo2')

tfbsGR<-GRanges(seqnames=BSwC$chr
                ,ranges=IRanges(start=BSwC$start,end=BSwC$end)
                ,strand='*'
                ,type=BSwC$type
                ,value=BSwC$value)

tfbsname<-unique(tfbsGR$type)


#######################################
#####   methylation infos
#######################################
prob_data<-fread('Data/Methylation/GSE62111_series_matrix.txt',header=T,sep = '\t',stringsAsFactors = F)

library(FDb.InfiniumMethylation.hg19)
# hm450 <- get450k()##
# save(hm450,file='hm450.Rdata')

# library(FDb.InfiniumMethylation.hg19)
# hm450 <- get450k()##
# save(hm450,file='hm450.Rdata')
load('Data/Methylation/hm450.Rdata')
y<-as.data.table(hm450)
y[,'ID_REF']<-names(hm450)
tmp<-prob_data[,c("ID_REF","GSM1519791")] #GM12878 methylation https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62111

names(tmp)<-c('ID_REF','v1')
prob_data_all<-left_join(tmp,y=y)

probGR<-GRanges(seqnames=(prob_data_all$seqnames)
                ,ranges=IRanges(start=prob_data_all$start,end=prob_data_all$start)
                ,strand='*'
                ,v1=prob_data_all$v1
                ,p1=0)
me <- probGR



####################################################################
#####          TSS information
####################################################################
tss<-read_tsv('Data/TSS/tss_RefGene.txt')
hks<-read_table2('Data/TSS/tss_HK_genes.txt',col_names = c('name','id'))
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
chr_info_1<-read_tsv('Data/Location/center.txt')
chr_info_2<-read_tsv('Data/Location/chr_size.txt',col_names = F)

chrnames<-c(paste('chr',1:22,sep=''),'chrX')

save(list=ls(),file='pre_info_08132020.Rdata')


