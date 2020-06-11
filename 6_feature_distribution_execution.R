#
#
#
#

gen_regions(test_ratio=3/10,BR=BR)
gen_data_1(region_GR)
save(anGR,file="anGR_compare_MCF7_MCF10A_T47D.Rdata")
load("anGR_compare_MCF7_MCF10A_T47D.Rdata")
#Large dataframe of features for TAD regions
#tmpsave <- anGR
anGR<-tmpsave
region_GR #contains the type of change and training and testing labeling

anGR <- as.data.frame(anGR)
regionGR <- as.data.frame(region_GR)

combine <- cbind(regionGR,anGR[,6:ncol(anGR)])
save(combine,file="anGR_compare_MCF7_MCF10A_T47D_combine_df.Rdata")




#####
#
#
# Start here
#
#
#
#####

load("anGR_compare_MCF7_MCF10A_T47D_combine_df.Rdata")

#anGR

features <- colnames(combine)[8:(ncol(combine)-1)] #remove chr location and to_cent
length(features)

combine[which(combine$change%in%c("MCF7only","T47Donly")),"change"] <- "canceronly"
classification <- c("MCF10Aonly","anyoverlap_all","canceronly")

train_data_1 <- combine[which(combine$train=="1"),]
test_data_1 <- combine[which(combine$train=="0"),]

x <- features
x <- c("H3K9ac_7","H3K9ac_10A","H3K4me1_7","H3K4me1_10A",
       "RAD21_7","RAD21_T47D","tfbs_RAD21",
       "tfbs_SMC3",
       "tfbs_CTCF","CTCF_7","CTCF_dmso_T47D","CTCF_10A",
       "tss_nc","tss_c","tss_hk","tss",
       "tfbs_STAT3","tfbs_CEBPB","tfbs_MAFK","tfbs_YY1",
       "H3K27me3_10A")
y <- "change"

train_data_1<-train_data_1[which(train_data_1$change%in%classification),c(y,x)]
test_data_1<-test_data_1[which(test_data_1$change%in%classification),c(y,x)]
train_data_1[,"change"]<-as.factor(train_data_1[,"change"])
test_data_1[,"change"]<-as.factor(test_data_1[,"change"])

#GBM multiclass http://uc-r.github.io/gbm_regression
library(h2o)
localH2O<-h2o.init(nthreads = 6,max_mem_size = "32G")
train.h2o <- as.h2o(train_data_1)
test.h2o <- as.h2o(test_data_1)
h2o.fit2<-h2o.gbm(x=x, y=y,training_frame = train.h2o,
                  nfolds = 5,
                  ntrees = 500,
                  max_depth = 10,
                  distribution = "AUTO")

h2o.fit2
perf <- h2o.performance(model=h2o.fit2,newdata=test.h2o)
perf
h2o.rmse(h2o.fit2, xval = TRUE)
vp <- h2o.varimp(h2o.fit2)
h2o.varimp_plot(h2o.fit2, num_of_features = 15)

save(h2o.fit2,vp,x, perf, file="anGR_compare_results_normal_vs_cancer_selectfeaturesmore.Rdata")
model_path <- h2o.saveModel(object=h2o.fit2, path=getwd(), force=TRUE)
print(model_path)

# jpeg(paste("normal_vs_breast_nocentro_ROC_curve.jpg",sep=''))
# plot(perf,type="roc")
# dev.off()

load("anGR_compare_results_normal_vs_cancer_nocentro.Rdata")
vp


tmp<-h2o.varimp(h2o.fit2)
var<-tmp[,1]
var <- sapply(strsplit(var,"_T47D"), `[`, 1)
var <- sapply(strsplit(var,"_7"), `[`, 1)# for MCF7
var <- sapply(strsplit(var,"_10A"), `[`, 1)# for MCF10a
var <- sapply(strsplit(var,"_dmso"), `[`, 1)# for MCF10a

tmp[,1]<-as.character(var)
uniquevar <- unique(var)
uniquevar[1:10]

vp_unique<-matrix(nrow=dim(tmp)[1],ncol=3)
rownames(vp_unique)<-uniquevar[1:dim(tmp)[1]]
colnames(vp_unique)<-c("relative_importance", "scaled_importance", "percentage")

for(i in 1:nrow(vp_unique)){
  vp_unique[i,"relative_importance"]<-sum(tmp[which(tmp[,1]==rownames(vp_unique)[i]),2])
  vp_unique[i,"scaled_importance"]<-sum(tmp[which(tmp[,1]==rownames(vp_unique)[i]),3])
  vp_unique[i,"percentage"]<-sum(tmp[which(tmp[,1]==rownames(vp_unique)[i]),4])
}

vp_unique[1:15,]
vp_unique_scaled <- sort(vp_unique[,"scaled_importance"],decreasing = T)
vp_unique_scaled <- vp_unique_scaled/vp_unique_scaled[1]
vp_forplot <- cbind(vp_unique[names(vp_unique_scaled)],vp_unique_scaled)

data<-as.data.frame(vp_forplot[1:15,])
data$'Feature'<-rownames(data)

# ggplot(data=data, aes(x=reorder(Feature,vp_unique_scaled),y=vp_unique_scaled)) +
#   geom_bar(position="dodge",stat="identity") + 
#   coord_flip() +
#   ylab("Scaled Importance") + 
#   xlab("Feature") +
#   ggtitle("Feature Importance Summary")


barplot(data$vp_unique_scaled,names.arg=data$Feature, horiz=T)
barplot(sort(data$vp_unique_scaled, decreasing = FALSE),names.arg=rev(data$Feature), horiz=T,
        xlab="Scaled Importance",
        main="Top 15 Feature Importance")
title(ylab="Features",mgp=c(7,1,0))
par(las=2)
