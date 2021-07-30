source("http://bioconductor.org/biocLite.R")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
biocLite("genefu")
biocLite("amap")
biocLite("breastCancerMAINZ")
biocLite("breastCancerTRANSBIG")

install.packages('R.utils')
install.packages('xtable')
install.packages('caret')
library(amap)
library(genefu)
library(xtable)
library(caret)
library(breastCancerMAINZ)
library(breastCancerTRANSBIG)
library(data.table)
library(R.utils)
setwd("E:/0_test/01_TNBC/0_data")

mRNA_exp_count <- fread("TCGA-BRCA_mRNA_counts.txt.gz",header = T,quote = "",sep="\t")
mRNA_exp_count <- as.data.frame(mRNA_exp_count)
miRNA_exp_count <- fread("TCGA-BRCA_counts.txt.gz",header = T,quote = "",sep="\t")
miRNA_exp_count <- as.data.frame(miRNA_exp_count)
sample_info <- fread("BRCA_survival.txt.gz",header = T,quote = "",sep="\t")
length(intersect(as.matrix(sample_info)[,1],colnames(mRNA_exp_count)))
sample_info_more <- fread("TCGA-BRCA_clinical.csv",header = T,quote = "",sep=",")
length(intersect(as.matrix(sample_info_more)[,1],substr(colnames(mRNA_exp_count),1,12)))

####### genefu pam50 #######
exp_mRNA_count <- mRNA_exp_count[,-1]
rownames(exp_mRNA_count) <- mRNA_exp_count[,1]
labels_sample <- as.numeric(substr(colnames(exp_mRNA_count),14,15))
samples_norm <- colnames(exp_mRNA_count)[which(labels_sample>=10)]
samples_dis <- colnames(exp_mRNA_count)[which(labels_sample<10)]
IDlist <- as.matrix(read.table("gene_ENSG2Symbol2ID.txt",header = F,quote = "",sep="\t"))
interID <- intersect(IDlist[,1],rownames(exp_mRNA_count))
rownames(IDlist) <- IDlist[,1]
EntrezGene.ID <- IDlist[interID,3]
ind_NA <- which(is.na(EntrezGene.ID))
interID_noNA <- interID[-ind_NA]
EntrezGene.ID_noNA <- EntrezGene.ID[-ind_NA]
ID_probe <- as.matrix(read.table("gene_probe.txt",header = F,quote = "",sep="\t",fill = T))
interID_probe <- intersect(interID_noNA,ID_probe[,4])
ID_probe_useful <- ID_probe[which(ID_probe[,4]%in%interID_probe),]
rownames(ID_probe_useful) <- ID_probe_useful[,4]
probes <- ID_probe_useful[interID_probe,1]
EntrezGene.ID_noNA_probe <- as.numeric(as.character(EntrezGene.ID_noNA[interID_probe]))
ddata <- t(exp_mRNA_count[interID_probe,samples_dis])
dannot <- data.frame(probe=probes,Gene.symbol=interID_probe,EntrezGene.ID=EntrezGene.ID_noNA_probe)
colnames(ddata) <- dannot$probe
PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=ddata,
                                annot=dannot,do.mapping=TRUE)
sample_class <- data.frame(samples=names(PAM50Preds$subtype),class=as.character(PAM50Preds$subtype))
sample_class <- rbind(as.matrix(sample_class),cbind(samples_norm,"norm"))
write.table(sample_class,"sample_class.txt",quote=F,row.names=F,col.names=T,sep = "\t")

setwd("E:\\0_test\\01_TNBC")
save.image(file = "allInfo.RData")
 load("allInfo.RData")
