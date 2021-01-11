# Get TAD boundaries for MCF7, MCF10A, and T47D
# Compare boundary for the three cell lines: gain, loss, any overlap, complete overlap
# Generate feature data for these regions
# Generate training and testing datasets and model alterations in TAD boundaries

library(plyr)
library(stringr)
library(readr)

load('Data/tads_bds.Rdata')

#MCF7 and MCF10A boundaries
bds<-ldply(f_b,function(x){x})
bds_mcf7<-bds[str_detect(bds$.id,'MCF7'),]
bds_mcf10a<-bds[str_detect(bds$.id,'MCF10a'),]

chrnames<-c(paste('chr',1:22,sep=''),'chrX')

bds_mcf7_chrkk<-bds_mcf7[str_detect(bds_mcf7$header,paste(chrnames[1],':',sep='')),]
chkk<-chrnames[1]
bds_mcf10a_chrkk<-bds_mcf10a[str_detect(bds_mcf10a$header,paste(chrnames[1],':',sep='')),]
chkk<-chrnames[1]
bdsGR_7<-GRanges(seqnames=chkk
                 ,ranges=IRanges(start=bds_mcf7_chrkk$start,end=bds_mcf7_chrkk$end)
                 ,strand='*')
bdsGR_10a<-GRanges(seqnames=chkk
                   ,ranges=IRanges(start=bds_mcf10a_chrkk$start,end=bds_mcf10a_chrkk$end)
                   ,strand='*')

for(kk in 2:23){
  bds_mcf7_chrkk<-bds_mcf7[str_detect(bds_mcf7$header,paste(chrnames[kk],':',sep='')),]
  chkk<-chrnames[kk]
  bds_mcf10a_chrkk<-bds_mcf10a[str_detect(bds_mcf10a$header,paste(chrnames[kk],':',sep='')),]
  chkk<-chrnames[kk]
  
  tmp1<-GRanges(seqnames=chkk
                 ,ranges=IRanges(start=bds_mcf7_chrkk$start,end=bds_mcf7_chrkk$end)
                 ,strand='*')
  bdsGR_7<-c(bdsGR_7, tmp1)
  tmp2<-GRanges(seqnames=chkk
                   ,ranges=IRanges(start=bds_mcf10a_chrkk$start,end=bds_mcf10a_chrkk$end)
                   ,strand='*')
  bdsGR_10a<-c(bdsGR_10a, tmp2)
}


#T47D boundaries
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


#GRanges
#bdsGR_7
#bdsGR_10a
#bdsGR_T47D


#Identify boundary regions that are gain or loss
#gain = in MCF7 but not in MCF10a
#loss = in MCF10a but not in MCf7

GR7<- bdsGR_7
GR10a <- bdsGR_10a
GRT47D <- bdsGR_T47D

#Overlap in MCF7 and T47D
BRoverlap_cancer <- GR7[findOverlaps(GR7, GRT47D,
                              type= "any",
                              select="all",
                              ignore.strand=FALSE)@from]
BRoverlap_cancer$type <- "anyoverlap_cancer"

BRoverlap_cancer100<- GR7[findOverlaps(GR7, GRT47D,
                                type= "equal",
                                select="all",
                                ignore.strand=FALSE)@from]
BRoverlap_cancer100$type<-"overlap_cancer"


#Overlap in all three cell lines
BRoverlap_all <- GR10a[findOverlaps(GR10a, BRoverlap_cancer,
                                     type= "any",
                                     select="all",
                                     ignore.strand=FALSE)@from]
BRoverlap_all$type <- "anyoverlap_all"

BRoverlap_all100<- GR10a[findOverlaps(GR10a, BRoverlap_cancer100,
                                       type= "equal",
                                       select="all",
                                       ignore.strand=FALSE)@from]
BRoverlap_all100$type<-"overlap_all"


#Comparing MCF7 and MCF10a

BRgain_710 <- GR7[!GR7 %over% GR10a,]
BRgain_710$type <- "gain_710"
BRloss_710 <- GR10a[!GR10a %over% GR7,]
BRloss_710$type <- "loss_710"

BRoverlap_710 <- GR7[findOverlaps(GR7, GR10a,
                              type= "any",
                              select="all",
                              ignore.strand=FALSE)@from]
BRoverlap_710$type <- "anyoverlap_710"

BRoverlap_710_100<- GR7[findOverlaps(GR7, GR10a,
                                type= "equal",
                                select="all",
                                ignore.strand=FALSE)@from]
BRoverlap_710_100$type<-"overlap_710"



#Comparing T47D and MCF10a

BRgain_T47D10 <- GRT47D[!GRT47D %over% GR10a,]
BRgain_T47D10$type <- "gain_T47D10"
BRloss_T47D10 <- GR10a[!GR10a %over% GRT47D,]
BRloss_T47D10$type <- "loss_T47D10"

BRoverlap_T47D10 <- GRT47D[findOverlaps(GRT47D, GR10a,
                              type= "any",
                              select="all",
                              ignore.strand=FALSE)@from]
BRoverlap_T47D10$type <- "anyoverlap_T47D10"

BRoverlap_T47D10_100<- GRT47D[findOverlaps(GRT47D, GR10a,
                                type= "equal",
                                select="all",
                                ignore.strand=FALSE)@from]
BRoverlap_T47D10_100$type<-"overlap_T47D10"



#MCF7, MCF10A, and T47D specific boundaries
BRgain_T47Donly <- GRT47D[!GRT47D %over% union(GR7,GR10a),]
BRgain_T47Donly$type <- "T47Donly"

BRgain_7only <- GR7[!GR7 %over% union(GRT47D,GR10a),]
BRgain_7only$type <- "MCF7only"

BRgain_10only <- GR10a[!GR10a %over% union(GRT47D,GR7),]
BRgain_10only$type <- "MCF10Aonly"


gr<-c(BRoverlap_cancer,BRoverlap_cancer100,BRoverlap_all,BRoverlap_all100,
      BRgain_710,BRloss_710,BRoverlap_710,BRoverlap_710_100,
      BRgain_T47D10,BRloss_T47D10,BRoverlap_T47D10,BRoverlap_T47D10_100,
      BRgain_T47Donly,BRgain_7only,BRgain_10only)

df <- data.frame(chr=seqnames(gr),
                 start=start(gr)-1,
                 end=end(gr)-1,
                 type=gr$type)

write.csv(df, file="br_gainloss_060820.csv", row.names=F,quote=F)

BR <- read.csv("br_gainloss_060820.csv",header=T,sep=',')
BR <- GRanges(seqnames=(BR$chr),
              ranges=IRanges(start=BR$start,end=BR$end),
              strand='*',
              change=BR$type)





# ##################################################################
# ##################################################################
# #
# # Preprocess feature information for MCF7 and MCF10A
# #
# ##################################################################
# ##################################################################
# 
# 
# setwd("Z:/SBMI_Houston/Jacqueline Chyr/TAD_alteration/R_codes_ZZG_for_TADs/0821")
# 
# ###################
# # ChIPseq data
# ####################
# 
# #H3K4me3_7, H4K20me1_7, H3K27ac_7, SMARCE1_7, FOXA1_7, CTCF_7, H3K27me3_7, H3K9me3_7, H3K4me2_7, SMARCA5_7, H3K4me1_7, H3K9ac_7
# bedf1<-read_table2('./mcf7_chip/List-of-Bed_narrowPeak-files(hg19).txt',col_names = T)
# bed_data1<-lapply(1:dim(bedf1)[1]
#                  ,function(i){read_table2(paste('./mcf7_chip/',bedf1$file_name[i],'.bed.gz',sep=''),col_names = F)}
# )
# 
# #EGR1_7, ELF1_7, EP300_7, GABPA_7, GATA3_7, H3K36me3_7, MAX_7, PML_7, POLR2A_7, RAD21_7, SIN3A_7, SRF_7
# bedf2<-read_table2('./mcf7_chip/List-of-Bam_files.txt',col_names = T)
# bed_data2<-lapply(1:dim(bedf2)[1]
#                   ,function(i){read_table2(paste('./mcf7_chip/bam_files/macs2/',bedf2$file_name[i],'_peaks.narrowPeak',sep=''),col_names = F)}
# )
# 
# #ChIPseq for MCF10A
# #broad H3K4me3_10A, narrow H3K27me3_10A, broad H3K79me2_10A, narrow CTCF_10A
# bedf3<-read_table2('./mcf10a_chip/List-of-Bed_narrowPeak-files(hg19).txt',col_names = T)
# bed_data3<-lapply(1:dim(bedf3)[1]
#                   ,function(i){read_table2(paste('./mcf10a_chip/',bedf3$file_name[i],sep=''),col_names = F)}
# )
# 
# #H3K4me1_10A, H3K4me3_10A, H3K9ac_10A, H3K9me3_10A, H3K27ac_10A, H3K27me3_10A, H3K36me3_10A, H3K79me2_10A
# bedf4<-read_table2('./mcf10a_chip/GSE85158_narrowPeak/List-of-Bed_narrowPeak-files(hg19).txt',col_names = T)
# bed_data4<-lapply(1:dim(bedf4)[1]
#                   ,function(i){read_table2(paste('./mcf10a_chip/GSE85158_narrowPeak/',bedf4$file_name[i],sep=''),col_names = F)}
# )
# 
# #ATACseq
# bedf5<-read_table2('./ATACseq/List-of-Bed_narrowPeak-files(hg19).txt',col_names = T)
# bed_data5<-lapply(1:dim(bedf5)[1]
#                   ,function(i){read_table2(paste('./ATACseq/',bedf5$file_name[i],sep=''),col_names = F)}
# )
# 
# names(bed_data1)<-bedf1$target #MCF7
# names(bed_data2)<-bedf2$target #MCF7
# names(bed_data3)<-bedf3$target #MCF10A broad and narrow (only CTCF is needed)
# names(bed_data4)<-bedf4$target #MCF10A
# names(bed_data5)<-bedf5$target #ATACseq
# 
# # peak_col_names<-c('chrom','start','end','name','score','strand','sig','pv','qv','peak','v1','v2','v3','v4')
# peak_col_names<-c('chrom','start','end','name','score','strand','sig','pv','qv','peak')
# 
# 
# 
# #######################################
# #####   methylation infos
# #######################################
# prob_data<-fread('./mcf7_chip/methylations/GSE71626_processed.txt',header=T,sep = '\t',stringsAsFactors = F)
# 
# # library(FDb.InfiniumMethylation.hg19)
# # hm450 <- get450k()##
# # save(hm450,file='hm450.Rdata')
# 
# # load('hm450.Rdata')
# # y<-as.data.table(hm450)
# # y[,'ID_REF']<-names(hm450)
# # tmp<-prob_data[,1:7]
# # names(tmp)<-c('ID_REF','v1','p1','v2','p2','v3','p3')
# # prob_data_all<-left_join(tmp,y=y)
# # 
# # probGR<-GRanges(seqnames=(prob_data_all$seqnames)
# #                 ,ranges=IRanges(start=prob_data_all$start,end=prob_data_all$start)
# #                 ,strand='*'
# #                 ,v1=prob_data_all$v1
# #                 ,p1=prob_data_all$p1
# #                 ,v2=prob_data_all$v2
# #                 ,p2=prob_data_all$p2
# #                 ,v3=prob_data_all$v3
# #                 ,p3=prob_data_all$p3)
# 
# load('me_data.Rdata')
# 
# 
# 
# 
# #################################
# #################################
# # Preprocess feature information for T47D
# #################################
# #################################
# 
# setwd("Z:/SBMI_Houston/Jacqueline Chyr/TAD_alteration/R_codes_T47D")
# 
# ###################
# # ChIPseq data
# ####################
# 
# bedf6<-read_table2('../Data/ChIPseq-out/T47D_hg19/List-of-narrowPeaks-T47D.txt',col_names = T)
# bed_data6<-lapply(1:dim(bedf6)[1]
#                  ,function(i){read_table2(paste('../Data/ChIPseq-out/T47D_hg19/',bedf6$file_name[i],'.bed.gz',sep=''),col_names = F)}
# )
# 
# 
# bedf7<-read_table2('../Data/ChIPseq-out/T47D_hg19/List-of-bam-T47D.txt',col_names = T)
# bed_data7<-lapply(1:dim(bedf7)[1]
#                   ,function(i){read_table2(paste('../Data/ChIPseq-out/T47D_hg19/bam_files/macs2/',bedf7$file_name[i],'_peaks.narrowPeak',sep=''),col_names = F)}
# )
# 
# 
# # bedf8<-read_table2('../Data/ChIPseq-out/T47D_hg38/List-of-narrowPeaks-hg38.txt',col_names = T)
# # chain<-import.chain('../Data/Annotations/liftOver/hg38ToHg19.over.chain')
# # for (i in 1:dim(bedf8)[1]){
# #   x1 <- read.table(paste("../Data/ChIPseq-out/T47D_hg38/",bedf8$file_name[i],".bed.gz",sep=''),header=F)
# #   x2 <- GRanges(seqnames=x1$V1,
# #                 ranges=IRanges(start=x1$V2,end=x1$V3),
# #                 strand='*',
# #                 V4=x1$V4,
# #                 V5=x1$V5,
# #                 V6=x1$V6,
# #                 V7=x1$V7,
# #                 V8=x1$V8)
# #   x3 <- liftOver(x2,chain)
# #   x4<-data.frame(unlist(x3),stringsAsFactors = F)
# #   write.table(x4,file=paste("../Data/ChIPseq-out/T47D_hg38/liftover-hg38-to-hg19/",bedf3$file_name[i],"_060920.hg19.bed",sep=''),quote=F,sep='\t',col.names=F, row.names=F)
# # }
# # rm(x1,x2,x3,x4)
# bed_data8<-lapply(1:dim(bedf8)[1]
#                   ,function(i){read_table2(paste('../Data/ChIPseq-out/T47D_hg38/liftover-hg38-to-hg19/',bedf8$file_name[i],'.hg19.bed',sep=''),col_names = F)})
# 
# peak_col_names<-c('chrom','start','end','name','score','strand','sig','pv','qv','peak')
# 
# names(bed_data6)<-bedf6$target #T47D
# names(bed_data7)<-bedf7$target #T47D
# names(bed_data8)<-bedf8$target #T47D
# 
# 
# #######################################
# #####   methylation infos
# #######################################
# #See gen_pre_info_01112020.R for T47D
# 
# load('me_data.Rdata')
# 
# 
# ###############################################
# ###     TFBS information (Transcription factors binding sites)
# ##############################################
# 
# setwd("Z:/SBMI_Houston/Jacqueline Chyr/TAD_alteration/R_codes_ZZG_for_TADs/0821")
# 
# BS<-fread('./mcf7_chip/BindingSites/wgEncodeRegTfbsClusteredV3.bed',sep = '\t',header = F)
# BSwC<-fread('./mcf7_chip/BindingSites/wgEncodeRegTfbsClusteredWithCellsV3.bed',sep = '\t',header = F)
# names(BSwC)<-c('chr','start','end','type','value','cell')
# names(BS)<-c('chr','start','end','type','value','NoC','cellNo1','cellNo2')
# 
# 
# tfbsGR<-GRanges(seqnames=BSwC$chr
#                 ,ranges=IRanges(start=BSwC$start,end=BSwC$end)
#                 ,strand='*'
#                 ,type=BSwC$type
#                 ,value=BSwC$value)
# 
# tfbsname<-unique(tfbsGR$type)
# 
# 
# ####################################################################
# #####          TSS information
# ####################################################################
# tss<-read_tsv('tss_RefGene.txt')
# hks<-read_table2('tss_HK_genes.txt',col_names = c('name','id'))
# tss<-tss[tss$chrom %in% chrnames,]
# 
# 
# #using the  unique row in tssfile 
# tmp<-duplicated(tss[,c('chrom','txStart','txEnd','name2')])
# tss<-tss[!tmp,]
# 
# #location of non-coding Tss
# loc_non_coding<-str_detect(tss$name2,'LINC*|LOC*|MIR*')
# 
# #location of coding Tss
# loc_coding<-!loc_non_coding
# 
# #location of HK Tss
# loc_hk<-tss$name2 %in% unique(hks$name)
# 
# ####deal with non_coding Tss with bins
# tssGR<-GRanges(seqnames=tss$chrom,
#                ranges=IRanges(start=tss$txStart,end=tss$txEnd),
#                strand=tss$strand)
# 
# ###################################################
# ###   location info 
# ##################################################
# ##get center and chr_size information for each chrs
# chr_info_1<-read_tsv('./mcf7_chip/center_chr/center.txt')
# chr_info_2<-read_tsv('./mcf7_chip/center_chr/chr_size.txt',col_names = F)
# 
# 
# 
# 
# 
# 
# setwd("Z:/SBMI_Houston/Jacqueline Chyr/TAD_alteration/Manuscript/BRCA-BR/After NAR review/After Journal of Genetics and Genomics review/Codes_for_manuscript_2020")
# 
# #save(list=ls(),file='pre_info_gain_loss_06082020.Rdata')


library(GenomicRanges)
library(data.table)

load('pre_info_gain_loss_06082020.Rdata')

BR <- read.csv("br_gainloss_060820.csv",header=T,sep=',')
BR <- GRanges(seqnames=(BR$chr),
              ranges=IRanges(start=BR$start,end=BR$end),
              strand='*',
              change=BR$type)

#6/10/2020 Plot distribution relative to TAD boundary center +/- 1 MB
region_GR<-bdsGR_T47D







#########################################
#########################################
#
# From source_funs_gainloss_0127_paper.r
# Generate training and testing dataset
#
#########################################
#########################################


# Randomly assign training and testing regions
chrnames<-c(paste('chr',1:22,sep=''),'chrX')

gen_regions<-function(test_ratio=3/10,BR=BR){
  region_GR<-BR
  
  # #10kb regions
  # final_range <- "final_ranges"
  # for(i in 1:length(region_GR)){
  #   end_range <- region_GR@ranges@start[i]+200000
  #   ranges <- end_range
  #   for(j in 1:20){
  #     end_range = end_range - 10000
  #     ranges <- c(end_range,ranges)
  #   }
  #   final_range<- c(final_range,ranges)
  # }
  # 
  # region_GR <- GRanges(seqnames=rep(as.character(region_GR@seqnames),each=21),
  #                 ranges=IRanges(start=as.numeric(final_range[-1]),end=as.numeric(final_range[-1])+10000),
  #                 strand='*',
  #                 change=rep(as.character(region_GR$change),each=21))
  # #end
  
  region_change<-unique(region_GR$change)
  for(i in 1:length(region_change)){
    assign(paste("region_",region_change[i],sep=''),region_GR[region_GR$change==region_change[i]])
  }
  
  # Randomize labeling
  #
  # region_bin_len = 200000
  # random_GR<-GRanges()
  # for (kk in 1:23){
  #   chrklen<-chr_info_2$X2[chr_info_2$X1%in% chrnames[kk]]
  #   chrkseq<-seq(from=1,by=region_bin_len,to=chrklen)
  #   bin_num<-length(chrkseq)
  #   tmpGR<-GRanges(seqnames = chrnames[kk]
  #                  ,ranges = IRanges(start=chrkseq[1:(bin_num-1)],end = chrkseq[2:bin_num])
  #                  ,strand = '*'
  #                  ,change = rep("noboundary",bin_num-1)
  #                  ,score = rep("0",bin_num-1)
  #   )
  #   random_GR<-append(random_GR,tmpGR)
  # }
  # '%!in%' <- function(x,y)!('%in%'(x,y))
  # region_random <- random_GR[sample(which(random_GR %!in% region_GR),size = 4107)]
  
  for(i in 1:length(region_change)){
    region_x <- eval(as.symbol(paste("region_",region_change[i],sep='')))
    samples<-length(region_x)
    tr <- rep(1,(samples))
    tr[sample(1:(samples),size = test_ratio*(samples))]<-0
    region_x$train<-tr
    assign(paste("region_",region_change[i],sep=''),region_x)
  }
  
  region_change_list<-paste("region_",region_change,sep='')
  tmp1<-sapply(region_change_list,as.symbol)
  tmp2<-sapply(tmp1,eval)
  region_GR <- c(tmp2[[1]],tmp2[[2]],tmp2[[3]],tmp2[[4]],tmp2[[5]],
                 tmp2[[6]],tmp2[[7]],tmp2[[8]],tmp2[[9]],tmp2[[10]],
                 tmp2[[11]],tmp2[[12]],tmp2[[13]],tmp2[[14]],tmp2[[15]])
  
  train_region <- subset(region_GR,region_GR$train==1)
  test_region<-subset(region_GR,region_GR$train!=1)
  return(list(region=region_GR,train=train_region,test=test_region))
}

#
#
#
#
#
#
#
#Generate feature values for those regions
#

gen_data_1<-function(region_GR){
  
  anGR <- region_GR
  bed_data1<-bed_data #change names
  bedf1<-bedfl #change names
  
  #Add ChIPseq data
  for(bd in c(1:2,4:8)){
    bed_data<-eval(as.symbol(paste("bed_data",bd,sep='')))
    for(i in 1:length(bed_data)){
      pt<-bed_data[[i]]
      names(pt)<-peak_col_names
      ptGR<-GRanges(seqnames=pt$chrom
                  ,ranges=IRanges(start=pt$start,end=pt$end)
                  ,strand='*'
                  ,sigv=pt$sig)
    
      tmp<-findOverlaps(ptGR,anGR,type='within',ignore.strand=T)
      tmp1<-data.table(pt=ptGR$sigv[tmp@from],id=tmp@to)
    
      mcols(anGR)[,names(bed_data)[[i]]]<-NA
      tmp2<-tmp1[,.(mv=mean(pt,na.rm=T)),by=id]
      mcols(anGR)[,names(bed_data)[[i]]][tmp2$id]<-tmp2$mv
    }
  }
  
  
  # for MCF10A ChIPseq. 
  # Length Class  Mode
  # H3K4me3_10A   9     tbl_df list
  # H3K27me3_10A 10     tbl_df list
  # H3K79me2_10A  9     tbl_df list
  # CTCF_10A     10     tbl_df list
  
  #ONLY CTCF
  for(i in 4){
    pt<-bed_data3[[i]]
    names(pt)<-peak_col_names
    ptGR<-GRanges(seqnames=pt$chrom
                  ,ranges=IRanges(start=pt$start,end=pt$end)
                  ,strand='*'
                  ,sigv=pt$sig)
    
    tmp<-findOverlaps(ptGR,anGR,type='within',ignore.strand=T)
    tmp1<-data.table(pt=ptGR$sigv[tmp@from],id=tmp@to)
    
    mcols(anGR)[,bedf3$target[i]]<-NA
    tmp2<-tmp1[,.(mv=mean(pt,na.rm=T)),by=id]
    mcols(anGR)[,bedf3$target[i]][tmp2$id]<-tmp2$mv
  }
  
  
  # tmp<-findOverlaps(probGR,anGR,type = 'within',ignore.strand=T)
  # prob_valid<-cbind(probGR$p1<0.05,probGR$p2<0.05,probGR$p3<0.05)
  # prob_v<-cbind(probGR$v1,probGR$v2,probGR$v3)
  # tmp1<-data.table(prob=(rowSums(prob_v*prob_valid)/rowSums(prob_valid))[tmp@from],id=tmp@to)
  # mcols(anGR)[,'me']<-NA
  # tmp2<-tmp1[,.(mv=mean(prob,na.rm=T)),by=id]
  # mcols(anGR)[,'me'][tmp2$id]<-tmp2$mv
  
  tmp<-findOverlaps(mcf10_me,anGR,type='within',ignore.strand=T)
  prob_valid <- mcf10_me$pvalue<0.05
  prob_v<-cbind(mcf10_me$beta)
  tmp1<-data.table(prob=(rowSums(prob_v*prob_valid)/rowSums(prob_valid))[tmp@from],id=tmp@to)
  mcols(anGR)[,'me_10A']<-NA
  tmp2<-tmp1[,.(mv=mean(prob,na.rm=T)),by=id]
  mcols(anGR)[,'me_10A'][tmp2$id]<-tmp2$mv
  
  tmp<-findOverlaps(mcf7_me,anGR,type='within',ignore.strand=T)
  prob_valid <- mcf7_me$pvalue<0.05
  prob_v<-cbind(mcf7_me$beta)
  tmp1<-data.table(prob=(rowSums(prob_v*prob_valid)/rowSums(prob_valid))[tmp@from],id=tmp@to)
  mcols(anGR)[,'me_7']<-NA
  tmp2<-tmp1[,.(mv=mean(prob,na.rm=T)),by=id]
  mcols(anGR)[,'me_7'][tmp2$id]<-tmp2$mv
  
  T47D_me$v1<-as.matrix(as.numeric(T47D_me$v1))
  T47D_me$p1<-as.matrix(as.numeric(T47D_me$p1))
  tmp<-findOverlaps(T47D_me,anGR,type='within',ignore.strand=T)
  prob_valid <- T47D_me$p1<0.05
  prob_v<-cbind(T47D_me$v1)
  tmp1<-data.table(prob=(rowSums(prob_v*prob_valid)/rowSums(prob_valid))[tmp@from],id=tmp@to)
  mcols(anGR)[,'me_T47D']<-NA
  tmp2<-tmp1[,.(mv=mean(prob,na.rm=T)),by=id]
  mcols(anGR)[,'me_T47D'][tmp2$id]<-tmp2$mv
  
  
  ###############################################
  ###     TFBS information (Transcription factors binding sites)
  ##############################################
  
  for(gene in tfbsname){
    tfbs_tmp<-tfbsGR[tfbsGR$type==gene,]
    tmp<-findOverlaps(tfbs_tmp,anGR,type = 'within',ignore.strand=T)
    tmp1<-data.table(v=tfbs_tmp$value[tmp@from],id=tmp@to)
    mcols(anGR)[,paste('tfbs',gene,sep='_')]<-NA
    tmp2<-tmp1[,.(mv=mean(v,na.rm=T)),by=id]
    mcols(anGR)[,paste('tfbs',gene,sep='_')][tmp2$id]<-tmp2$mv
  }
  
  ####################################################################
  #####          TSS information
  ####################################################################
  
  #non coding genes
  tmpGR<-GRanges(seqnames=tss$chrom[loc_non_coding],
                 ranges=IRanges(start=tss$txStart[loc_non_coding],end=tss$txEnd[loc_non_coding]),
                 strand=tss$strand[loc_non_coding])
  tmp<-findOverlaps(tmpGR,anGR,type='any',ignore.strand=T)
  mcols(anGR)[,'tss_nc']<-NA
  tmp1<-data.table(from=tmp@from,to=tmp@to)
  tmp2<-tmp1[,n:=length(from),by=to]
  
  mcols(anGR)[,'tss_nc'][tmp2$to]<-tmp2$n

  ####deal with coding Tss with bins
  tmpGR<-GRanges(seqnames=tss$chrom[loc_coding],
                 ranges=IRanges(start=tss$txStart[loc_coding],end=tss$txEnd[loc_coding]),
                 strand=tss$strand[loc_coding])
  tmp<-findOverlaps(tmpGR,anGR,type='any',ignore.strand=T)
  mcols(anGR)[,'tss_c']<-NA
  tmp1<-data.table(from=tmp@from,to=tmp@to)
  tmp2<-tmp1[,n:=length(from),by=to]
  mcols(anGR)[,'tss_c'][tmp2$to]<-tmp2$n
  
  ####deal with hk Tss with bins
  tmpGR<-GRanges(seqnames=tss$chrom[loc_hk],
                 ranges=IRanges(start=tss$txStart[loc_hk],end=tss$txEnd[loc_hk]),
                 strand=tss$strand[loc_hk])
  tmp<-findOverlaps(tmpGR,anGR,type='any',ignore.strand=T)
  mcols(anGR)[,'tss_hk']<-NA
  tmp1<-data.table(from=tmp@from,to=tmp@to)
  tmp2<-tmp1[,n:=length(from),by=to]
  mcols(anGR)[,'tss_hk'][tmp2$to]<-tmp2$n
  
  #### all TSS
  tmp<-findOverlaps(tssGR,anGR,type='any',ignore.strand=T)
  mcols(anGR)[,'tss']<-NA
  tmp1<-data.table(from=tmp@from,to=tmp@to)
  tmp2<-tmp1[,n:=length(from),by=to]
  mcols(anGR)[,'tss'][tmp2$to]<-tmp2$n
  
  
  ###################################################
  ###   location info 
  ##################################################
  ##get center and chr_size information for each chrs
  ## 0 = close to the centromere, 1 = further away
  i=1
  chr_end<-c()
  chr_cent<-c()
  for(i in 1:length(anGR)){
    chr_end[i]<-chr_info_2$X2[chr_info_2$X1==as.character(anGR@seqnames[i])]
    chr_cent[i]<-chr_info_1$chromStart[chr_info_1$chrom==as.character(anGR@seqnames[i])&chr_info_1$type=='centromere']
  }
  
  to_cent_tmp<-(chr_cent-anGR@ranges@start+anGR@ranges@width/2)
  anGR$to_cent<-to_cent_tmp/ifelse(to_cent_tmp>0,chr_cent,chr_cent-chr_end)
  # anGR$to_tele<-(anGR@ranges@start+anGR@ranges@width/2)/chr_end
  
  
  
  #
  #
  #
  #
  #
  #
  #
  #
  # Differences in signals
  #
  #
  
  
  # anGR$CTCF_diff <- anGR$CTCF_7 - anGR$CTCF_10A
  # anGR$H3K4me1_diff <- anGR$H3K4me1_7 - anGR$H3K4me1_10A
  # anGR$H3K4me3_diff <- anGR$H3K4me3_7 - anGR$H3K4me3_10A
  # anGR$H3K9ac_diff <- anGR$H3K9ac_7 - anGR$H3K9ac_10A
  # anGR$H3K27ac_diff <- anGR$H3K27ac_7 - anGR$H3K27ac_10A
  # anGR$H3K27me3_diff <- anGR$H3K27me3_7 - anGR$H3K27me3_10A
  # anGR$H3K36me3_diff <- anGR$H3K36me3_7 - anGR$H3K36me3_10A
  # 
  # anGR$ATACseq_diff <- anGR$ATACseq_7 - anGR$ATACseq_10A
  # 
  # anGR$me_diff <- anGR$me_7 - anGR$me_10A
  
  return(anGR)
}
save(anGR,file="anGR_compare_MCF7_MCF10A_T47D_200kb.Rdata")

library(gplots)
#boxplots
df <- as.data.frame(anGR)
feature <- "tfbs_SMC3"
a<-df[which(df$change=="MCF10Aonly"),"CTCF_10A"]
b<-df[which(df$change=="MCF7only"),"CTCF_7"]
c<-df[which(df$change=="T47Donly"),"CTCF_dmso_T47D"]

a<-df[which(df$change=="MCF10Aonly"),feature]
b<-df[which(df$change%in%c("T47Donly","MCF7only")),feature]
c<-df[which(df$change=="overlap_all"),feature]

boxplot(a,b,c,
        names=c("Normal","Cancer","Conserved"),
        main = "SMC3 TFBS",
        xlab=("TAD Boundaries in Cell lines"),
        ylab="Number of Binding Sites")
par(las=1)
#TADBoundariesDistribution_CTCFChIPseq


