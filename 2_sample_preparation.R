# These are codes to prepare samples for PredTAD
# Last edited August 2020
# feature_data_preparation, previously known as gen_pre_info
# Data is saved in all_data.Rdata


load('pre_info_08132020.Rdata') #preprocessed ChIPseq and epigenomic and genomic data


region_bin_len=10000 #width of each genomic region
test_ratio=3/10 #ratio of the testing set
neighbors=10 #number of neghboring bins to consider

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

###############
###############
##
## Extract the data for these samples
##
###############
###############

bdsGR <-bdsGR_GM12878 #select the cell line you want

gen_data_1<-function(train_region=region_GR,by_len=10000,n_neighbors=10){
  ana_list<-list()
  for(kk in 1:23){
    
    #construct data for modeling
    chkk<-chrnames[kk]
    bdsGR_chrkk<-bdsGR[which(bdsGR@seqnames == chrnames[kk])]
    ana_data<-data.frame(start=seq(from=0,to=max(bdsGR@ranges@start+bdsGR@ranges@width),by=by_len)+1
                         ,end=seq(from=by_len,to=max(bdsGR@ranges@start+bdsGR@ranges@width)+by_len,by=by_len))
    
    #using GRanges ****************************************************************
    anGR<-GRanges(seqnames=chkk
                  ,ranges=IRanges(start=ana_data$start,end=ana_data$end)
                  ,strand='*')
    
    tmp<-findOverlaps(anGR,train_region,type='within',ignore.strand=T)
    anGR<-anGR[tmp@from]
    
    ##################################################
    ###     add chipSeq informations 
    ##################################################
    
    #for histone and TF ChIPseq data
    for(i in 1:length(bed_data1)){
      pt<-bed_data1[[i]]
      names(pt)<-peak_col_names
      ptGR<-GRanges(seqnames=pt$chrom
                    ,ranges=IRanges(start=pt$start,end=pt$end)
                    ,strand='*'
                    ,sigv=pt$sig)
      
      tmp<-findOverlaps(ptGR,anGR,type='within',ignore.strand=T)
      tmp1<-data.table(pt=ptGR$sigv[tmp@from],id=tmp@to)
      
      mcols(anGR)[,bedfl$target[i]]<-NA
      tmp2<-tmp1[,.(mv=mean(pt,na.rm=T)),by=id]
      mcols(anGR)[,bedfl$target[i]][tmp2$id]<-tmp2$mv
    }
    
    #normalize average chIPseq peaks
    for(i in 1:length(bed_data1)){
      x<-anGR@elementMetadata[,i]
      iqr <- IQR(x,na.rm=TRUE)
      max <- iqr +quantile(x, probs=c(0.75),na.rm=TRUE)
      anGR@elementMetadata[,i] <- anGR@elementMetadata[,i]/max
      anGR@elementMetadata[which(anGR@elementMetadata[,i] > 1),i] <- 1
    }
    
    #################################################
    ###     add methylation values
    #################################################
    tmp<-findOverlaps(me,anGR,type='within',ignore.strand=T)
    prob_valid <- me$p1<0.05 #p-value
    prob_v<-cbind(me$v1) #methylation beta value
    tmp1<-data.table(prob=(rowSums(prob_v*prob_valid)/prob_valid)[tmp@from],id=tmp@to)
    mcols(anGR)[,'me']<-NA
    tmp2<-tmp1[,.(mv=mean(prob,na.rm=T)),by=id]
    mcols(anGR)[,'me'][tmp2$id]<-tmp2$mv
    
    
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
    
    tmpGR<-GRanges(seqnames=tss$chrom[loc_non_coding],
                   ranges=IRanges(start=tss$txStart[loc_non_coding],end=tss$txEnd[loc_non_coding]),
                   strand=tss$strand[loc_non_coding])
    
    #bins
    tmp<-findOverlaps(tmpGR,anGR,type='any',ignore.strand=T)
    # tmp
    
    mcols(anGR)[,'tss_nc']<-NA
    tmp1<-data.table(from=tmp@from,to=tmp@to)
    tmp2<-tmp1[,n:=length(from),by=to]
    mcols(anGR)[,'tss_nc'][tmp2$to]<-tmp2$n
    
    ####deal with coding Tss with bins
    tmpGR<-GRanges(seqnames=tss$chrom[loc_coding],
                   ranges=IRanges(start=tss$txStart[loc_coding],end=tss$txEnd[loc_coding]),
                   strand=tss$strand[loc_coding])
    
    #bins
    tmp<-findOverlaps(tmpGR,anGR,type='any',ignore.strand=T) #### tmpGR? or tssGR?
    # tmp
    
    mcols(anGR)[,'tss_c']<-NA
    tmp1<-data.table(from=tmp@from,to=tmp@to)
    tmp2<-tmp1[,n:=length(from),by=to]
    mcols(anGR)[,'tss_c'][tmp2$to]<-tmp2$n
    
    ####deal with hk Tss with bins
    tmpGR<-GRanges(seqnames=tss$chrom[loc_hk],
                   ranges=IRanges(start=tss$txStart[loc_hk],end=tss$txEnd[loc_hk]),
                   strand=tss$strand[loc_hk])
    
    #bins
    tmp<-findOverlaps(tmpGR,anGR,type='any',ignore.strand=T)
    # tmp
    
    mcols(anGR)[,'tss_hk']<-NA
    tmp1<-data.table(from=tmp@from,to=tmp@to)
    tmp2<-tmp1[,n:=length(from),by=to]
    mcols(anGR)[,'tss_hk'][tmp2$to]<-tmp2$n
    
    
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
    
    
    #####################################################
    ### set the boundaries to be 1
    tmp<-findOverlaps(anGR,bdsGR,type='any',ignore.strand=T)
    mcols(anGR)$bdsGR<-0
    mcols(anGR)$bdsGR[tmp@from]<-1
    
    ana_list[[kk]]<-anGR
    # ana_list[[kk]]<-anGR[c(-1:-10,-length(anGR):(-length(anGR)+9))]
    cat('kk=',kk,'\n')
    
  }
  
  
  ######################################
  ##  gen neighbours 
  vars<-names(mcols(ana_list[[1]]))
  vars<-setdiff(vars,c('bdsGR'))
  
  tictoc::tic()
  
  ana_dt_old<-as.data.table(ldply(ana_list,
                                  function(x){cbind(chr=as.vector(x@seqnames),
                                                    as.data.frame(x@ranges)[,1:2],
                                                    as.data.frame(mcols(x)))}))
  setkey(ana_dt_old,chr,start,end)
  ana_dt_old$id<-1:dim(ana_dt_old)[1]
  ana_dt_old[,paste(rep(vars,each=n_neighbors),'L',rep(1:n_neighbors,length(vars)),sep='_'):=rep(0,dim(ana_dt_old)[1])]
  ana_dt_old[,paste(rep(vars,each=n_neighbors),'R',rep(1:n_neighbors,length(vars)),sep='_'):=rep(0,dim(ana_dt_old)[1])]
  
  rs<-dim(ana_dt_old)[1]
  for(j in 1:n_neighbors){
    for(i in 1:length(vars)){
      ana_dt_old[[paste(vars[i],'L',j,sep='_')]]<-c(rep(NA,j),ana_dt_old[[vars[i]]][1:(rs-j)])
      ana_dt_old[[paste(vars[i],'R',j,sep='_')]]<-c(ana_dt_old[[vars[i]]][(j+1):rs],rep(NA,j))
    }
  }
  tictoc::toc()
  
  ana_dt_old$del<-F
  del_bin<-function(del,id){
    n<-length(id)
    del[c(1:n_neighbors,(n-n_neighbors+1):n)]<-T
    return(del)
  }
  ana_dt_old[,.(del=del_bin(del,id)),by=chr]->t1
  ana_dt_old$del<-t1$del
  ana_dt_old_bak<-ana_dt_old
  ana_dt_old<-ana_dt_old[del!=T]
  ana_dt_old$del<-NULL
  
  ana_dt_old[is.na(ana_dt_old)]<-0
  ana_dt_old$bds<-factor(ana_dt_old$bds,levels = c(0,1),labels = c(0,1))
  ana_dt_old$chr<-factor(ana_dt_old$chr,levels=chrnames,labels=chrnames)
  return(ana_dt_old)
}


all_data_1 <- gen_data_1(region_GR,by_len=region_bin_len,n_neighbors=neighbors)
save(all_data_1,region_GR,file="all_data_1_GM12878.Rdata")


