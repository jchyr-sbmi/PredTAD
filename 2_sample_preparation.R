# These are codes to prepare samples for PredTAD
# Last edited May 26, 2020
# feature_data_preparation, previously known as gen_pre_info
# Data is saved in all_data.Rdata


load('pre_info_12102019.Rdata') #preprocessed ChIPseq and epigenomic and genomic data


region_bin_len=10000 #width of each genomic region
test_ratio=3/10 #ratio of the testing set
n_neighbors=10 #number of neghboring bins to consider

region_GR<-GRanges()
for (kk in 1:23){
  chrklen<-chr_info_2$X2[chr_info_2$X1%in% chrnames[kk]]
  chrkseq<-seq(from=1,by=region_bin_len,to=chrklen)
  bin_num<-length(chrkseq)
  tr<-rep(1,(bin_num-1))
  tr[sample(1:(bin_num-1),size = test_ratio*(bin_num-1))]<-0
  tmpGR<-GRanges(seqnames = chrnames[kk]
                 ,ranges = IRanges(start=chrkseq[1:(bin_num-1)],end = chrkseq[2:bin_num])
                 ,strand = '*'
                 ,train= tr
  )
  region_GR<-append(region_GR,tmpGR)
}

chr_info_tc <- chr_info_1[which(chr_info_1$type %in% c("centromere","telomere")),] # identify telomere and centromere regions
chr_info_tc <- GRanges(seqnames = chr_info_tc$chrom,
                       ranges = IRanges(start=chr_info_tc$chromStart, end = chr_info_tc$chromEnd),
                       strand = "*",
                       type = chr_info_tc$type)

tmp_ov <- findOverlaps(region_GR, chr_info_tc,type="any")
region_GR$telocentro <- "0"
region_GR$telocentro[tmp_ov@from]<- "1"

all_data_1 <- gen_data_1(region_GR,by_len=region_bin_len,k_neighbour=n_neighbors)
save(all_data_1,region_GR,file="all_data_1_01092020.Rdata")








