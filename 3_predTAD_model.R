rm(list=ls())
gc()
library(tidyverse)
library(stringr)
library(data.table)
library(ROSE)
library(randomForest)
library(ROCR)
library(h2o)
library(readr)
library(minfi)
library(GenomicRanges)
library(plyr)
library(tictoc)
library(h2o)
localH2O<-h2o.init(nthreads = 6,max_mem_size = "32G")
h2o.removeAll()

# PredTAD codes
# Run feature_data_preparation.R first to generate preprocessed data.


# all_data_1<-read.table("Data/all_data.txt",sep='\t',header=T,stringsAsFactors = F)
# sample_regions <-read.table("Data/sample.txt",sep='\t',header=T,stringsAsFactors = F)
# 
# df<-all_data_1
# df[,'bds_7']<-factor(df[,'bds_7'])
# df[,'bds_10a']<-factor(df[,'bds_10a'])
# df[,'bds_T47D']<-factor(df[,'bds_T47D'])
# df[,'chr']<-factor(df[,'chr'])
# save(df,file="Z:/SBMI_Houston/Jacqueline Chyr/TAD_alteration/Manuscript/BRCA-BR/After NAR review/After Journal of Genetics and Genomics review/Codes_for_manuscript_2020/Data/all_data.Rdata")

load("Data/all_data.Rdata")

# test_ratio = .30
# tr<-rep(1,nrow(all_data_1))
# tr[sample(1:nrow(all_data_1),size = test_ratio*nrow(all_data_1)-1)]<-0
# 
# all_data_1 <- cbind(tr,all_data_1)
# save(tr,file="tr.Rdata")
load("Data/tr.Rdata")

all_data_1<-cbind(tr,df)
# all_data_1<- as.data.frame(all_data_1)
# all_data_1[,5]<-as.character(all_data_1[,5])
# all_data_1[,6]<-as.character(all_data_1[,6])
# 
# sample_regions[,2]<-as.character(sample_regions[,2])
# sample_regions[,3]<-as.character(sample_regions[,3])

#colnames(sample_regions)[1]<-"chr"
#write.table(sample_regions,file="Data/sample.txt",sep='\t',col.names=T,row.names=F,quote=F)

#df<-merge(sample_regions,all_data_1,by.1="chr",by.2="start",by.3="end",all=FALSE,no.dups = TRUE)

#write.table(all_data_1,file="Z:/SBMI_Houston/Jacqueline Chyr/TAD_alteration/Manuscript/BRCA-BR/After NAR review/After Journal of Genetics and Genomics review/Codes_for_manuscript_2020/Data/all_data.txt",quote=F,col.names=T,row.names=F,sep='\t')
'%!in%' <- function(x,y)!('%in%'(x,y))

chr_odd <- c("chr1","chr3","chr5","chr7","chr9","chr11","chr13","chr15","chr17","chr19","chr21")
chr_even <-c("chr2","chr4","chr6","chr8","chr10","chr12","chr14","chr16","chr18","chr20","chr22")

train_data_1<- all_data_1[which(all_data_1$chr%in%chr_odd),]
test_data_1 <- all_data_1[which(all_data_1$chr%in%chr_even),]

train_data_1<- all_data_1[which(all_data_1$chr%!in%c("chr1","chr8","chr19")),]
test_data_1 <- all_data_1[which(all_data_1$chr%in%c("chr1","chr8","chr19")),]

train_data_1<- all_data_1[which(all_data_1$tr==1),] #70%
test_data_1 <- all_data_1[which(all_data_1$tr==0),] #30%

dim(train_data_1)
dim(test_data_1)

data_h2o<-as.h2o(train_data_1)
test_h2o<-as.h2o(test_data_1)

# train_data_1<- all_data_1[which(region_GR$train==1 & region_GR$telocentro=="0"),]
# test_data_1 <- all_data_1[which(region_GR$train==0 & region_GR$telocentro=="0"),]
# dim(train_data_1)
# dim(test_data_1)


#Select your features
tfbs <- colnames(train_data_1)[str_detect(colnames(train_data_1),'tfbs')]
tss <- colnames(train_data_1)[str_detect(colnames(train_data_1),'tss')]
to_cent <- colnames(train_data_1)[str_detect(colnames(train_data_1),'to_cent')]
me10a <- colnames(train_data_1)[str_detect(colnames(train_data_1),'me_10A')]
me7 <- colnames(train_data_1)[str_detect(colnames(train_data_1),'me_7')]
atacseq10a <- colnames(train_data_1)[str_detect(colnames(train_data_1),'ATACseq_10A')]
atacseq7 <- colnames(train_data_1)[str_detect(colnames(train_data_1),'ATACseq_7')]
T47D <- colnames(train_data_1)[str_detect(colnames(train_data_1),'T47D')]
T47D <- T47D[-(which(T47D == "bds_T47D"))]
mcf10a <- colnames(train_data_1)[str_detect(colnames(train_data_1),'_10A')]
mcf10a <- mcf10a[mcf10a %!in% c(me10a,atacseq10a)]
mcf7 <- colnames(train_data_1)[str_detect(colnames(train_data_1),'_7')]
mcf7 <- mcf7[mcf7 %!in% c(tfbs,tss,to_cent,me10a,me7,atacseq10a,atacseq7,mcf10a)]
mcf7 <- mcf7[-(which(mcf7 == "bds_7"))]
  

#Select the cell line and the features
y <- "bds_T47D"
x <- c('chr',tfbs,tss,to_cent,T47D)

y <- "bds_7"
x <- c('chr',tfbs,tss,to_cent,me7,atacseq7,mcf7)

y <- "bds_10a"
x <- c('chr',tfbs,tss,to_cent,me10a,atacseq10a,mcf10a)

top15_MCF10A <- c("chr",
           colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_SMC3")==TRUE)], #2, 21
           colnames(train_data_1)[which(colnames(train_data_1)=="tfbs_CTCF")],#3, 41
           colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_CTCF_")==TRUE)], #3, 41
           colnames(train_data_1)[which(startsWith(colnames(train_data_1),"H3K9ac_10A")==TRUE)], #4
           colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_RAD21")==TRUE)], #5
           colnames(train_data_1)[which(startsWith(colnames(train_data_1),"H3K27me3_10A")==TRUE)], #6
           colnames(train_data_1)[which(startsWith(colnames(train_data_1),"H3K9me3_10A")==TRUE)], #7
           colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tss_nc")==TRUE)], #8
           colnames(train_data_1)[which(startsWith(colnames(train_data_1),"to_cent")==TRUE)], #9
           colnames(train_data_1)[which(startsWith(colnames(train_data_1),"H3K36me3_10A")==TRUE)], #10
           colnames(train_data_1)[which(startsWith(colnames(train_data_1),"me_10A")==TRUE)], #11
           colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_ZNF143")==TRUE)], #12
           colnames(train_data_1)[which(colnames(train_data_1)=="tfbs_FOS")],
           colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_FOS_")==TRUE)], #13
           colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_EP300")==TRUE)], #14
           colnames(train_data_1)[which(startsWith(colnames(train_data_1),"H3K4me1_10A")==TRUE)] #15
           )

top15_MCF7 <- c("chr",
                  colnames(train_data_1)[which(startsWith(colnames(train_data_1),"CTCF_7")==TRUE)], #2, 21
                  colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tss_nc")==TRUE)],
                  colnames(train_data_1)[which(colnames(train_data_1)=="tfbs_CTCF")],#3, 41
                  colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_CTCF_")==TRUE)], #4
                  colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_SMC3")==TRUE)], #5
                  colnames(train_data_1)[which(startsWith(colnames(train_data_1),"RAD21_7")==TRUE)], #6
                  colnames(train_data_1)[which(startsWith(colnames(train_data_1),"to_cent")==TRUE)], #7
                  colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_RAD21")==TRUE)], #8
                  colnames(train_data_1)[which(startsWith(colnames(train_data_1),"ATACseq_7")==TRUE)], #9
                  colnames(train_data_1)[which(startsWith(colnames(train_data_1),"H3K9me3_7")==TRUE)], #10
                  colnames(train_data_1)[which(startsWith(colnames(train_data_1),"me_7")==TRUE)], #11
                  colnames(train_data_1)[which(startsWith(colnames(train_data_1),"H3K4me1_7")==TRUE)], #12
                  colnames(train_data_1)[which(startsWith(colnames(train_data_1),"POLR2A_7")==TRUE)], #13
                  colnames(train_data_1)[which(startsWith(colnames(train_data_1),"GATA3_7")==TRUE)], #14
                  colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_CEBPB")==TRUE)] #15
                )

top15_T47D <- c("chr",
                colnames(train_data_1)[which(startsWith(colnames(train_data_1),"CTCF_T47D")==TRUE)], #2, 21
                colnames(train_data_1)[which(startsWith(colnames(train_data_1),"to_cent")==TRUE)], #3, 41
                colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_RAD21")==TRUE)], #4
                colnames(train_data_1)[which(colnames(train_data_1)=="tfbs_CTCF")],#3, 41
                colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_CTCF_")==TRUE)], #5
                colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_SMC3")==TRUE)], #6
                colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tss_c")==TRUE)], #7
                colnames(train_data_1)[which(startsWith(colnames(train_data_1),"me_T47D")==TRUE)], #8
                colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_EP300")==TRUE)], #9
                colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_CEBPB")==TRUE)], #10
                colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_FOS")==TRUE)], #11
                colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_MAFK")==TRUE)], #12
                colnames(train_data_1)[which(startsWith(colnames(train_data_1),"ATACseq_T47D")==TRUE)], #13
                colnames(train_data_1)[which(startsWith(colnames(train_data_1),"tfbs_POLR2A")==TRUE)], #14
                colnames(train_data_1)[which(startsWith(colnames(train_data_1),"GATA3_T47D")==TRUE)] #15
                )







data_h2o<-as.h2o(train_data_1)
test_h2o<-as.h2o(test_data_1)

#1,8,19
#GBM model 
y<-"bds_7"
x<-c('chr',tfbs,tss,to_cent,me7,atacseq7,mcf7)

gbm <- h2o.gbm(x=x,y=y,training_frame = data_h2o, balance_classes = T, ntrees = 500, max_depth = 10,nfolds = 5,distribution="bernoulli")
vp <- h2o.varimp(gbm)
perf <- h2o.performance(model=gbm,newdata=test_h2o)
h2o.auc(gbm)

stats <- gbm@model$training_metrics@metrics$max_criteria_and_metric_scores
save(gbm,vp, perf, stats, file="m_compare_chr1_8_19_MCF7_training.Rdata")
model_path <- h2o.saveModel(object=gbm, path=getwd(), force=TRUE)
print(model_path)

pred.GBM<-h2o.performance(model=gbm,newdata=test_h2o)
jpeg(paste("chr1_8_19_MCF7_ROC_curve.jpg",sep=''))
plot(pred.GBM,type="roc")
dev.off()

test_data_1 <- all_data_1[which(all_data_1$chr%in%c("chr1")),]
#test_data_1 <- all_data_1[which(all_data_1$chr%in%c("chr8")),]
#test_data_1 <- all_data_1[which(all_data_1$chr%in%c("chr19")),]

test_h2o<-as.h2o(test_data_1)
perf <- h2o.performance(model=gbm,newdata=test_h2o)
save(perf, file="m_compare_chr1_8_19_MCF7_testing_chr1.Rdata")

#testing
test_data_1 <- all_data_1[which(all_data_1$chr%in%c("chr1")),]
test_data_1 <- all_data_1[which(all_data_1$chr%in%c("chr8")),]
test_data_1 <- all_data_1[which(all_data_1$chr%in%c("chr21")),]

test_h2o<-as.h2o(test_data_1)
perf <- h2o.performance(model=gbm,newdata=test_h2o)
save(perf, file="m_compare_chr1_8_19_MCF7_testing_chr21.Rdata")
perf




#rename saved model




#GBM model 
y<-"bds_10a"
x<-top15_MCF10A

data_h2o<-as.h2o(train_data_1[,c(y,x)])
test_h2o<-as.h2o(test_data_1[,c(y,x)])

gbm <- h2o.gbm(x=x,y=y,training_frame = data_h2o, balance_classes = T, ntrees = 500, max_depth = 10,nfolds = 5,distribution="bernoulli")
vp <- h2o.varimp(gbm)
perf <- h2o.performance(model=gbm,newdata=test_h2o)
h2o.auc(gbm)

stats <- gbm@model$training_metrics@metrics$max_criteria_and_metric_scores
save(gbm,vp, perf, stats, file="m_compare_top15_MCF10A_training.Rdata")
model_path <- h2o.saveModel(object=gbm, path=getwd(), force=TRUE)
print(model_path)

pred.GBM<-h2o.performance(model=gbm,newdata=test_h2o)
jpeg(paste("Top15_MCF10A_ROC_curve.jpg",sep=''))
plot(pred.GBM,type="roc")
dev.off()

#rename saved model





##############################6/6/20

#GBM model 
y<-"bds_7"
x<-top15_MCF7[1:190] #10
x<-top15_MCF7[1:85] #5

#data_h2o<-as.h2o(train_data_1[,c(y,x)])
#test_h2o<-as.h2o(test_data_1[,c(y,x)])

gbm <- h2o.gbm(x=x,y=y,training_frame = data_h2o, balance_classes = T, ntrees = 500, max_depth = 10,nfolds = 5,distribution="bernoulli")
vp <- h2o.varimp(gbm)
perf <- h2o.performance(model=gbm,newdata=test_h2o)
h2o.auc(gbm)

stats <- gbm@model$training_metrics@metrics$max_criteria_and_metric_scores
save(gbm,vp, perf, stats, file="m_compare_top5_withchr_MCF7_training.Rdata")
model_path <- h2o.saveModel(object=gbm, path=getwd(), force=TRUE)
print(model_path)
print("top5")

pred.GBM<-h2o.performance(model=gbm,newdata=test_h2o)
jpeg(paste("Top5_withchr_MCF7_ROC_curve.jpg",sep=''))
plot(pred.GBM,type="roc")
dev.off()

tmp<-vp
var<-tmp[,1] #variable names
var <- sapply(strsplit(var,"_L_"), `[`, 1)
var <- sapply(strsplit(var,"_R_"), `[`, 1)
var <- sapply(strsplit(var,"_7"), `[`, 1)

tmp[,1]<-as.character(var)
uniquevar <- unique(var)
uniquevar[1:10]
length(uniquevar)

vp_unique<-matrix(nrow=length(uniquevar),ncol=3)
rownames(vp_unique)<-uniquevar[1:length(uniquevar)]
colnames(vp_unique)<-c("relative_importance", "scaled_importance", "percentage")

for(i in 1:nrow(vp_unique)){
  vp_unique[i,"relative_importance"]<-sum(tmp[which(tmp[,1]==rownames(vp_unique)[i]),2])
  vp_unique[i,"scaled_importance"]<-sum(tmp[which(tmp[,1]==rownames(vp_unique)[i]),3])
  vp_unique[i,"percentage"]<-sum(tmp[which(tmp[,1]==rownames(vp_unique)[i]),4])
}

vp_unique[1:length(uniquevar),]
vp_unique_scaled <- sort(vp_unique[,"scaled_importance"],decreasing = T)
vp_unique_scaled <- vp_unique_scaled/vp_unique_scaled[1]
vp_forplot <- cbind(vp_unique[names(vp_unique_scaled)],vp_unique_scaled)

data<-as.data.frame(vp_forplot[1:length(uniquevar),])
data$'Feature'<-rownames(data)
data$Feature <- factor(data$Feature, levels = data$Feature[order(data$vp_unique_scaled)])

library(ggplot2)
ggplot(data=data, aes(x=Feature,y=vp_unique_scaled)) +
  geom_bar(position="dodge",stat="identity") + 
  coord_flip() +
  ylab("Scaled Importance") + 
  ggtitle("Feature Importance Summary for MCF7")








#GBM model 
y<-"bds_T47D"
x<-top15_T47D

gbm <- h2o.gbm(x=x,y=y,training_frame = data_h2o, balance_classes = T, ntrees = 500, max_depth = 10,nfolds = 5,distribution="bernoulli")
vp <- h2o.varimp(gbm)
perf <- h2o.performance(model=gbm,newdata=test_h2o)
h2o.auc(gbm)

stats <- gbm@model$training_metrics@metrics$max_criteria_and_metric_scores
save(gbm,vp, perf, stats, file="m_compare_top15_T47D_training.Rdata")
model_path <- h2o.saveModel(object=gbm, path=getwd(), force=TRUE)
print(model_path)

pred.GBM<-h2o.performance(model=gbm,newdata=test_h2o)
jpeg(paste("Top15_T47D_ROC_curve.jpg",sep=''))
plot(pred.GBM,type="roc")
dev.off()

