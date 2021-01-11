## 
## Analyze boundaries gained and lost between MCF10A and MCF7:
##
## 1. See if any estrogen-related genes or oncogenes or BR survival-related genes are 
##    located near the altered boundaries.
## 2. Check if any of those genes have predicted enhancers across the boundary
## 3. See if there's any changes in histone modification / ChIP-seq peaks
##
##
rm(list=ls());gc()
library(readr)
library(GenomicRanges)

## file for MCF7 and MCF10A boundaries
BR <- read.csv("BR-mcf7_10.csv",header=T,sep=',')
BR$chi <- paste("chr", BR$chi, sep="")

## locations where boundary is gained or lost
BR.gain <- BR[which(BR$change=='gain'),]
BR.loss <- BR[which(BR$change=='loss'),]

BR.gain.loc<-GRanges(seqnames = BR.gain$chi, IRanges(start=BR.gain$start,end=BR.gain$end),strand = '*')#,name=geneloc$name2)
BR.loss.loc<-GRanges(seqnames = BR.loss$chi, IRanges(start=BR.loss$start,end=BR.loss$end),strand = '*')#,name=geneloc$name2)

## List of oncogenes or tumor suppressors: see gene_expression_code_new.R
## http://ncg.kcl.ac.uk/cancer_genes.php
## http://www.uniprot.org/uniprot/?query=keyword:oncogene&fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score

oncogenes <- read.table("Gene-Expression/DEG.Cancerlistall.txt")

## List of estrogen related genes
## From codes_for_MCF7_E2/EBG_E2.R 
## MCF7-E2-12hr-GSE11352.txt logFC > 1 and adj-p < 0.2
## MCF7-17BE2-2hr-GSE55922-NURSA.txt logFC > 1 and p-val < 0.05

load("../../codes_for_MCF7_E2/MCF7_E2_genelist.RData")
E2_genes <- genelist

## List of differentiated genes (cell line) MCF7 vs MCF10A
## List of differentiated genes (TCGA) 
## Luminal A patients: ER+, PR+/-, HER2-
## 315 lumA samples vs 40 matched normal samples

cell.exp <- read_tsv("Gene-Expression/GSE71862_DeSeq2_MCF7_vs_MCF10A_foldchange.txt",col_names = T)
Targetgenes.cell <- unique(cell.exp[abs(cell.exp$`MCF7/MCF10A log2FoldChange`)>1 & 
                                      cell.exp$pvalue<0.05,"genename"])
DEG.cell.genes <- Targetgenes.cell$genename


#RET expression barplot
cell.exp.raw <- read_tsv("GSE71862_MCF7_MCF10A_RSEM_expectedcounts.txt",col_names=T)
RET <- cell.exp.raw[which(cell.exp.raw$gene=="RET"),]








load("Gene-Expression/Topgenes.TCGA.05.order.Rdata")
Topgenes.TCGA.05.order.no0 <- Topgenes.TCGA.05.order[Topgenes.TCGA.05.order[,"lum.med"]!="0" & 
                                                       Topgenes.TCGA.05.order[,"matlum.med"]!="0" & 
                                                       Topgenes.TCGA.05.order[,"matnor.med"]!="0" &
                                                       Topgenes.TCGA.05.order[,"allno.med"]!="0",]
# dim(Topgenes.TCGA.05.order)
# 4545   21
# dim(Topgenes.TCGA.05.order.no0)
# 3947   21
DEG.TCGA.genes <- sapply(str_split(rownames(Topgenes.TCGA.05.order.no0),'\\|'),function(x){x[1]})


length(E2_genes)
length(DEG.cell.genes)
length(DEG.TCGA.genes)



## Gene location
geneloc<-read_tsv("Gene-Expression/tss_RefGene.txt",col_names=T)
geneloc<-GRanges(seqnames = geneloc$chrom,ranges = IRanges(start=geneloc$txStart,end=geneloc$txEnd),strand = '*',name=geneloc$name2)
#geneloc<-geneloc[order(geneloc$name,-geneloc$width),]
geneloc<-geneloc[!duplicated(geneloc$name),]
#View(geneloc)





query <- geneloc
query <- geneloc[geneloc$name%in%DEG.cell.genes]
query <- geneloc[geneloc$name%in%E2_genes]

ext <- 1e+06
BR.gain.loc.ext<-GRanges(seqnames = BR.gain$chi, IRanges(start=BR.gain$start-ext,end=BR.gain$end+ext),strand = '*')
BR.loss.loc.ext<-GRanges(seqnames = BR.loss$chi, IRanges(start=BR.loss$start-ext,end=BR.loss$end+ext),strand = '*')

subject <- BR.gain.loc.ext
subject <- BR.loss.loc.ext


ov<-findOverlaps(query,subject,type='any',ignore.strand=T)
ov.genes<-query[unique(ov@from),]
length(ov.genes)


#For gene ontology analysis
write.table(ov.genes$name,file="genes_BRgain_1MB.txt",sep=',',quote=F,col.name=F,row.name=F)
write.table(ov.genes$name,file="genes_BRloss_1MB_DEGcelline.txt",sep=',',quote=F,col.name=F,row.name=F)
