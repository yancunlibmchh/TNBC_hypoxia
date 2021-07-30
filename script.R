library(RColorBrewer)

rm(list = ls())
options(stringsAsFactors = F)
source('Z://projects/codes/mg_base.R')
plotKMCox <- function(dat){
  colnames(dat)=c('time','status','groups')
  sdf<-survdiff(Surv(time,status) ~ groups,data=dat)
  print((sdf))
  p<-pchisq(sdf$chisq,length(sdf$n)-1,lower.tail=FALSE)
  sf<-survfit(Surv(time,status) ~ groups,data=dat)
  colKm=rainbow(length(sf$strata))
  plot(sf, mark.time = TRUE,col=colKm,xlab=paste("Survival time in day","\np=",round(p,5)),ylab = "Survival probabilities",main="Method Kaplan Meier")
  legend('topright',paste0(gsub('groups=','',names(sf$strata)),'(N=',sdf$n,')'), col = colKm,
         lty = c(1,1, 1, 1),lwd=c(1,1,1,1),merge = TRUE,cex = 0.8)
  return(p)
}
plotKMCox_1=function(dat,n){
  library(survival)
  library(ggsci)
  mypal = pal_jama(alpha = 0.7)(7)
  colnames(dat)=c('time','status','groups')
  sdf<-survdiff(Surv(time,status) ~ groups,data=dat)
  #print((sdf))
  p<-pchisq(sdf$chisq,length(sdf$n)-1,lower.tail=FALSE)
  sf<-survfit(Surv(time,status) ~ groups,data=dat)
  colKm=c(mypal[2],mypal[1],mypal[3],mypal[4])
  plot(sf, mark.time = TRUE,col=colKm,xlab=paste("Survival time in day","\np=",round(p,5)),ylab = "Survival probabilities",main=n)
  legend('topright',paste0(gsub('groups=','',names(sf$strata)),'(N=',sdf$n,')'), col = colKm,
         lty = c(1,1, 1, 1),lwd=c(1,1,1,1),merge = TRUE,cex = 0.8)
  return(p)
}


library(stringr)
library(dplyr)
GSE104193_cli <- read.delim('origin_datas/GSE104193/GSE104193_series_matrix.txt',header = T, stringsAsFactors = F, check.names = F)
GSE104193_cli <- GSE104193_cli[, c(2,10,11,13)]
GSE104193_cli <- GSE104193_cli[GSE104193_cli$Sample_characteristics_ch1 == 'cell line: MDA-MB-231', ]
GSE104193_cli$type <- str_split_fixed(GSE104193_cli$Sample_characteristics_ch1.1, ":", 2)[,2]
GSE104193_cli$code <- str_split_fixed(GSE104193_cli$Sample_characteristics_ch1.2, ":", 2)[,2]
GSE104193_cli <- GSE104193_cli[, -c(3,4)]
GSE104193_cli <- GSE104193_cli[GSE104193_cli$type == " Hypoxia" | GSE104193_cli$type == " Normoxia", ]
GSE104193_cli <- arrange(GSE104193_cli, type)
GSE104193_cli$type <-  trimws(GSE104193_cli$type, which = c("both"), whitespace = "[ \t\r\n]")
GSE104193_cli$code <-  trimws(GSE104193_cli$code, which = c("both"), whitespace = "[ \t\r\n]")

GSE104193_exp <- read.delim('origin_datas/GSE104193/GSE104193_exp_mRNA.txt',header = T, stringsAsFactors = F, check.names = F, row.names = 1)
GSE104193_exp <- GSE104193_exp[, GSE104193_cli$code]
colnames(GSE104193_exp) <- GSE104193_cli$Sample_geo_accession
GSE104193_exp <- apply(GSE104193_exp, 2, function(x) {round(x, 0)})



library(DESeq2)
eset=GSE104193_exp[, GSE104193_cli$Sample_geo_accession]
group_list <- c(rep('Hypoxia', 4),rep('Normoxia', 4))
colData <- data.frame(row.names=colnames(eset), group_list)

dds <- DESeqDataSetFromMatrix(eset, colData, design= ~ group_list)
keep<- rowMeans(counts(dds))>=5
dds.keep<-dds[keep,]
dds.keep<-DESeq(dds.keep)
res= results(dds.keep)

al=cbind(row.names(res), res$log2FoldChange, res$pvalue, res$padj)
colnames(al) <- c('genes', 'log2FC', 'pvalue', 'FDR')
class(al)
al <- crbind2DataFrame(al)
al <- na.omit(al)
GSE104193_diff <- al[abs(al$log2FC) > log2(1.5) & al$FDR < 0.05, ]
GSE104193_diff$log2FC <- -GSE104193_diff$log2FC
write.table(GSE104193_diff, file = 'files/GSE104193_diff_genes.txt',
            sep = '\t', row.names = F, quote = F)

GSE104193_diff_genes <- GSE104193_diff$genes
GSE104193_diff_up_genes <- GSE104193_diff[GSE104193_diff$log2FC > 0, ]$genes
GSE104193_diff_dn_genes <- GSE104193_diff[GSE104193_diff$log2FC < 0, ]$genes


GSE104193.deg.dif <- which(abs(as.numeric(al$log2FC)) > log2(1.5) & as.numeric(al$FDR) < 0.05)
pdf('PDFs/GSE104193_diff_gene_volcano.pdf',width = 5,height = 5)
dif <- al[GSE104193.deg.dif,]
t1 <- al
log2FC <- -as.numeric(t1$log2FC)
fdr <- -log2(as.numeric(t1$FDR))
col <- ifelse(fdr> -log2(0.05),ifelse(log2FC > log2(1.5),'#DC143C',ifelse(log2FC < -log2(1.5),'#4169E1','grey')),'grey')
plot(log2FC,fdr,pch=16,col=col,
     ylim=c(0,60),xlim=c(-4, 4),
     xlab='log2FC',
     ylab='-log2(FDR)',
     main=paste0('Hypoxia / Normoxia'))    
legend("topleft",box.col="black",
       c(paste0("Up :",as.numeric(table(col)[2])),
         paste0("Down :",as.numeric(table(col)[1]))),
       col=c("#DC143C","#4169E1"),pch=16,cex=1)
dev.off()


library(pheatmap)

dat <- eset[c(GSE104193_diff_up_genes, GSE104193_diff_dn_genes), 
            GSE104193_cli$Sample_geo_accession]

annotation_col  <- data.frame(Sample = factor(GSE104193_cli$type))

rownames(annotation_col)=colnames(dat)

bk=unique(c(seq(-1.5,1.5, length=100)))
pdf('PDFs/GSE104193_diff_gene_heatmap.pdf',width = 5,height = 5)
pheatmap::pheatmap(dat, scale = 'row', 
         breaks = bk,
         annotation_col = annotation_col,
         cluster_cols = F, cluster_rows = F,
         show_rownames = F, show_colnames = F,
         gaps_col = 4,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()



GSE33950_cli <- read.delim('origin_datas/GSE33950_family/SampleInfo_contrib1-GPL570.txt',header = T, stringsAsFactors = F, check.names = F)
GSE33950_cli <- GSE33950_cli[1:8, c("Accession", "Title")]
GSE33950_cli$type <- c(rep('Normoxia', 4),
                       rep('Hypoxia', 4))
GSE33950_cli <- arrange(GSE33950_cli, type)

GSE33950_exp <- read.delim('origin_datas/GSE33950_family/GSE33950_exp.txt',header = T, stringsAsFactors = F, check.names = F, row.names = 1)
boxplot(GSE33950_exp)
GSE33950_exp <- GSE33950_exp[, GSE33950_cli$Accession]
boxplot(GSE33950_exp)



library(limma)
package.version("limma")
group_list <- factor(GSE33950_cli$type)
design <- model.matrix(~group_list)
fit <- lmFit(GSE33950_exp, design)
fit <- eBayes(fit)
options(digits = 4)
topTable(fit,coef=2,adjust='BH')
deg <- topTable(fit,coef=2,adjust='BH',number = Inf)
deg$genes <- rownames(deg)

GSE33950_diff <- deg[abs(deg$logFC) > log2(1.5) & deg$adj.P.Val < 0.05, ]
GSE33950_diff$logFC <- -GSE33950_diff$logFC
write.table(GSE33950_diff, file = 'files/GSE33950_diff_genes.txt',
            sep = '\t', row.names = F, quote = F)


GSE33950_diff_genes <- GSE33950_diff$genes


GSE33950_diff_up_genes <- GSE33950_diff[GSE33950_diff$logFC > 0, ]$genes
GSE33950_diff_dn_genes <- GSE33950_diff[GSE33950_diff$logFC < 0, ]$genes


GSE33950.deg.dif <- which(abs(as.numeric(deg$logFC)) > log2(1.5) & as.numeric(deg$adj.P.Val) < 0.05)
pdf('PDFs/GSE33950_diff_gene_volcano.pdf',width = 5,height = 5)
dif <- deg[GSE33950.deg.dif,]
t1 <- deg
log2FC <- -as.numeric(t1$logFC)
fdr <- -log2(as.numeric(t1$adj.P.Val))
col <- ifelse(fdr> -log2(0.05),ifelse(log2FC > log2(1.5),'#DC143C',ifelse(log2FC < -log2(1.5),'#4169E1','grey')),'grey')
plot(log2FC,fdr,pch=16,col=col,
     ylim=c(0,30),xlim=c(-3, 3),
     xlab='log2FC',
     ylab='-log2(FDR)',
     main=paste0('Hypoxia / Normoxia'))    
legend("topright",box.col="black",
       c(paste0("Up :",as.numeric(table(col)[2])),
         paste0("Down :",as.numeric(table(col)[1]))),
       col=c("#DC143C","#4169E1"),pch=16,cex=1)
dev.off()



library(pheatmap)

dat <- GSE33950_exp[c(GSE33950_diff_up_genes, GSE33950_diff_dn_genes), 
            GSE33950_cli$Accession]

annotation_col  <- data.frame(Sample = factor(GSE33950_cli$type))

rownames(annotation_col)=colnames(dat)

bk=unique(c(seq(-1.5,1.5, length=100)))
pdf('PDFs/GSE33950_diff_gene_heatmap.pdf',width = 5,height = 5)
pheatmap::pheatmap(dat, scale = 'row', 
                   breaks = bk,
                   annotation_col = annotation_col,
                   cluster_cols = F, cluster_rows = F,
                   show_rownames = F, show_colnames = F,
                   gaps_col = 4,
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()



hypoxia_genes <- read.csv('origin_datas/S1_Table.hypoxia-genes.txt', header = F)
hypoxia_genes <- hypoxia_genes$V1





library(ggplot2)
library(VennDiagram)
venn2 <- venn.diagram(list(HRGs = hypoxia_genes, 
                           GSE33950 = GSE33950_diff_genes,
                           GSE104193 = GSE104193_diff_genes),
                      filename=NULL,
                      fill = c("#619CFF", "#00BA38", "#F8766D"))    
grid.draw(venn2)
dev.off()
ggsave(plot = venn2,
       filename = 'PDFs/dHRGsVenn.pdf',
       width = 5, height = 5, device = cairo_pdf)


GSE33950_dHRGs <- intersect(hypoxia_genes, GSE33950_diff_genes)
GSE104193_dHRGs <- intersect(hypoxia_genes, GSE104193_diff_genes)

dHRGs <- unique(c(GSE33950_dHRGs, GSE104193_dHRGs))

dHRGs_res <- enrichmentORA(dHRGs,
                           mp_dbs=c('pathway_KEGG',
                                    'geneontology_Biological_Process',
                                    'geneontology_Cellular_Component',
                                    'geneontology_Molecular_Function'))
dHRGs_res_filter <- dHRGs_res[dHRGs_res$FDR < 0.05, ]

pdf('PDFs/dHRGsEnrichment.pdf', width = 18, height = 8)
dotplot_batch(dHRGs_res_filter, 
              dbs =c('geneontology_Biological_Process',
                     'geneontology_Cellular_Component',
                     'geneontology_Molecular_Function',
                     'pathway_KEGG'),
              top=10, FDR = T)
dev.off()

write.table(dHRGs_res, file = 'files/dHRGsEnrichment.txt', 
            sep = '\t', row.names = F, quote = F)
table(dHRGs_res$DB)


GSE103091_cli <- GSE103091$Sample
GSE103091_exp <- GSE103091$Exp$GPL570_54675_Data_col2
GSE103091_anno <- GSE103091$Anno$GPL570

GSE103091_cli <- crbind2DataFrame(GSE103091_cli)
GSE103091_exp <- crbind2DataFrame(GSE103091_exp)
GSE103091_anno <- crbind2DataFrame(GSE103091_anno)

table(GSE103091_cli$Desc)

GSE103091_cli <- GSE103091_cli[, c('Acc', "os (days)", "death")]
GSE103091_cli$`os (days)`
GSE103091_cli <- GSE103091_cli[GSE103091_cli$`os (days)` != 'NULL', ]
GSE103091_cli <- crbind2DataFrame(GSE103091_cli)
colnames(GSE103091_cli) <- c('Samples', 'OS.time', 'OS')
rownames(GSE103091_cli) <- GSE103091_cli$Samples


rownames(GSE103091_anno) <- GSE103091_anno$V1
length(intersect(rownames(GSE103091_exp), GSE103091_anno$V1))

GSE103091_exp$genes <- GSE103091_anno[rownames(GSE103091_exp), ]$V11

GSE103091_exp <- aggregate(.~genes, data=GSE103091_exp, mean) 
GSE103091_exp <- GSE103091_exp[-c(1), ]
rownames(GSE103091_exp) <- GSE103091_exp$genes
GSE103091_exp <- GSE103091_exp[, -c(1)]

write.table(GSE103091_exp,
            file = 'origin_datas/GSE103091_exp.txt',
            sep = '\t', quote = F)

GSE103091_cli
write.table(GSE103091_cli,
            file = 'origin_datas/GSE103091_cli.txt',
            sep = '\t', quote = F)

intersect(dHRGs, rownames(GSE103091_exp))
setdiff(dHRGs, rownames(GSE103091_exp))
intersect('AK3', rownames(GSE103091_exp))

GSE103091_dHRGs <- intersect(dHRGs, rownames(GSE103091_exp))

intersect(GSE103091_dHRGs, rownames(GSE103091_exp))


GSE103091_exp_log2 <- t(log2(GSE103091_exp[, GSE103091_cli$Samples] + 1))
boxplot(t(GSE103091_exp_log2))

GSE103091_sig_cox_data <- cbind(GSE103091_cli, GSE103091_exp_log2)
GSE103091_sig_cox_data <- GSE103091_sig_cox_data[, -c(1)]

cox_univariant_gene_regr<-function(myd,colums){
  library(survival)
  Surv(as.numeric(myd$OS.time),as.numeric(myd$OS))->myd.surv
  c()->univar_anova_p;
  c()->univar_coxph_HR;
  c()->univar_coxph_low95;
  c()->univar_coxph_high95;
  c()->univar_coxph_logtest;
  c()->univar_coxph_sctest;
  c()->univar_coxph_waldtest;
  c()->fpkm.mean;
  c()->fpkm.median;
  colnames(myd)[colums]->myd.names;
  for(i in myd.names){
    mean(myd[,i],na.rm=T)->tmp.mean;
    median(myd[,i],na.rm=T)->tmp.median;
    c(fpkm.mean,tmp.mean)->fpkm.mean;
    c(fpkm.median,tmp.median)->fpkm.median;
    as.formula(paste("myd.surv~",i))->tmp.formula;
    coxph(formula=tmp.formula,data=myd)->tmp.coxph;
    summary(tmp.coxph)->tmp.coxph.summary;
    c(univar_anova_p,tmp.coxph.summary$coefficients[,5])->univar_anova_p;
    c(univar_coxph_HR,tmp.coxph.summary$coefficients[,2])->univar_coxph_HR;
    c(univar_coxph_low95,tmp.coxph.summary$conf.int[,3])->univar_coxph_low95;
    c(univar_coxph_high95,tmp.coxph.summary$conf.int[,4])->univar_coxph_high95;
    c(univar_coxph_logtest,tmp.coxph.summary$logtest[3])->univar_coxph_logtest;
    c(univar_coxph_sctest,tmp.coxph.summary$sctest[3])->univar_coxph_sctest;
    c(univar_coxph_waldtest,tmp.coxph.summary$waldtest[3])->univar_coxph_waldtest;	
  }
  data.frame("gName"=myd.names,
             "Pvalue"=univar_anova_p,
             "HR"=univar_coxph_HR,
             "Low(95%CI)"=univar_coxph_low95,
             "High(95%CI)"=univar_coxph_high95,
             "Logrank"=univar_coxph_logtest,
             "Sctest"=univar_coxph_sctest,
             "Waldtest"=univar_coxph_waldtest,
             "fpkm_median"=fpkm.median,
             "fpkm_mean"=fpkm.mean)->myd.coxph.df;
  myd.coxph.df[!is.na(myd.coxph.df$Logrank),]->myd.coxph.df;
  myd.coxph.df[order(myd.coxph.df$Logrank),]->myd.coxph.df;
  return(myd.coxph.df);
}
draw_forestplot<-function(fp.coxph){
  fp.coxph[,c("gName","Pvalue","HR","Low.95.CI.","High.95.CI.")]->fp.table;
  fp.table[order(fp.table$HR,decreasing=F),]->fp.table
  ggplot(data=fp.table,
         aes(x=gName,y=HR,ymin=Low.95.CI.,
             ymax=High.95.CI.)
  )->p;
  p+geom_point(aes(size=-log(Pvalue)),
               color=brewer.pal(11,"Spectral")[4])+geom_hline(yintercept=1,lty=2)+coord_flip()+xlab("Shared genes")+ylab("95% Confidence Interval")+geom_errorbar(aes(ymin=Low.95.CI.,ymax=High.95.CI.),width=0.5,cex=0.5,lty=2)+theme_bw()->p;
  p+scale_x_discrete(limits=fp.table$gName)->p;
  return(p);
}

GSE103091_genes_index <- which(colnames(GSE103091_sig_cox_data)%in%GSE103091_dHRGs)

GEO_exp_symbol_processed.coxph <- cox_univariant_gene_regr(GSE103091_sig_cox_data,GSE103091_genes_index)
GSE103091.forestplot.coxph <- GEO_exp_symbol_processed.coxph
GSE103091.forestplot.coxph <- GSE103091.forestplot.coxph[GSE103091.forestplot.coxph$Pvalue<0.05,]
dim(GSE103091.forestplot.coxph)
GSE103091.forestplot <- draw_forestplot(GSE103091.forestplot.coxph)
GSE103091.forestplot

GSE103091.fp.coxph <- GEO_exp_symbol_processed.coxph
GSE103091.fp.table <- GSE103091.fp.coxph[,c("gName","Pvalue","HR","Low.95.CI.","High.95.CI.")]
GSE103091.fp.table <- GSE103091.fp.table[order(GSE103091.fp.table$HR,decreasing=F),]
genes_order <- as.character(GSE103091.fp.table$gName)
rownames(GSE103091.fp.coxph) <- GSE103091.fp.coxph[,1]
genes_order_CI <- GSE103091.fp.coxph[genes_order[length(genes_order):1],2:5]
genes_order_CI_sig <- genes_order_CI[which(genes_order_CI[,1]<0.05),]
write.table(genes_order_CI,"files/GSE103091_genes_order_CI.txt",quote=F,row.names=T,col.names=T,sep="\t")
write.table(genes_order_CI_sig,"files/GSE103091_genes_order_CI_p0.05.txt",quote=F,row.names=T,col.names=T,sep="\t")




library(data.table)
cbioportal_tnbc_samples <- read.csv('origin_datas/cbioportal/TNBC_samples.txt', 
                                    header = T, sep = '\t', stringsAsFactors = F)
cbioportal_tnbc_samples <- cbioportal_tnbc_samples$PATIENT_ID


cbioportal_tnbc_cli <- read.delim('origin_datas/cbioportal/data_clinical_patient_OS.txt', header = T,stringsAsFactors = F, row.names = 1)
cbioportal_tnbc_cli <- na.omit(cbioportal_tnbc_cli)
colnames(cbioportal_tnbc_cli) <- c('OS.time', 'OS')
cbioportal_tnbc_cli$OS.time <- cbioportal_tnbc_cli$OS.time * 30
cbioportal_tnbc_cli$OS <- ifelse(cbioportal_tnbc_cli$OS == 'LIVING', 0, 1)




library(data.table)
cbioportal_exp <- fread('origin_datas/cbioportal/data_expression_median.txt',
                        sep = '\t',stringsAsFactors = F, data.table = F)
rownames(cbioportal_exp) <- cbioportal_exp$Hugo_Symbol
cbioportal_exp <- cbioportal_exp[, -c(1,2)]

intersect(rownames(cbioportal_exp), dHRGs)
setdiff(dHRGs, rownames(cbioportal_exp))


cbioportal_exp_dHRGs <- intersect(dHRGs, rownames(cbioportal_exp))
intersect(rownames(cbioportal_exp), cbioportal_exp_dHRGs)

cbioportal_tnbc_samples <- intersect(cbioportal_tnbc_samples,
                                     intersect(rownames(cbioportal_tnbc_cli),
                                               colnames(cbioportal_exp)))


cbioportal_tnbc_cli <- cbioportal_tnbc_cli[cbioportal_tnbc_samples, ]
head(cbioportal_tnbc_cli)

cbioportal_exp_filter <- cbioportal_exp[cbioportal_exp_dHRGs, cbioportal_tnbc_samples]
boxplot(cbioportal_exp_filter)

cbioportal_exp_filter <- t(cbioportal_exp_filter)
cbioportal_exp_filter <- cbind(cbioportal_tnbc_cli, cbioportal_exp_filter)
cbioportal_exp_filter <- cbioportal_exp_filter[cbioportal_exp_filter$OS.time > 0, ]



library(ggplot2)
cbioportal_genes_index <- which(colnames(cbioportal_exp_filter)%in%cbioportal_exp_dHRGs)

cbioportal_exp_symbol_processed.coxph <- cox_univariant_gene_regr(cbioportal_exp_filter,cbioportal_genes_index)
cbioportal.forestplot.coxph <- cbioportal_exp_symbol_processed.coxph
cbioportal.forestplot.coxph <- cbioportal.forestplot.coxph[cbioportal.forestplot.coxph$Pvalue<0.05,]
dim(cbioportal.forestplot.coxph)
cbioportal.forestplot <- draw_forestplot(cbioportal.forestplot.coxph)
cbioportal.forestplot

cbioportal.fp.coxph <- cbioportal_exp_symbol_processed.coxph
cbioportal.fp.table <- cbioportal.fp.coxph[,c("gName","Pvalue","HR","Low.95.CI.","High.95.CI.")]
cbioportal.fp.table <- cbioportal.fp.table[order(cbioportal.fp.table$HR,decreasing=F),]
genes_order <- as.character(cbioportal.fp.table$gName)
rownames(cbioportal.fp.coxph) <- cbioportal.fp.coxph[,1]
genes_order_CI <- cbioportal.fp.coxph[genes_order[length(genes_order):1],2:5]
genes_order_CI_sig <- genes_order_CI[which(genes_order_CI[,1]<0.05),]
write.table(genes_order_CI,"files/cbioportal_genes_order_CI.txt",quote=F,row.names=T,col.names=T,sep="\t")
write.table(genes_order_CI_sig,"files/cbioportal_genes_order_CI_p0.05.txt",quote=F,row.names=T,col.names=T,sep="\t")


library(data.table)
TCGA_tnbc <- fread('Z:/TCGA/Matrix/mRNA_TPM_Symbol/Merge_TCGA-BRCA_TPM.txt',sep = '\t',
                   stringsAsFactors = F, data.table = F)
rownames(TCGA_tnbc) <- TCGA_tnbc$Tag
TCGA_tnbc <- TCGA_tnbc[-c(1), ]

TCGA_tnbc_cli <- read.delim('origin_datas/TCGA/TCGA_BRCA_sample.txt', header = T,
                            stringsAsFactors = F)
TCGA_tnbc_cli <- TCGA_tnbc_cli[TCGA_tnbc_cli$breast_carcinoma_estrogen_receptor_status == 'Negative' & TCGA_tnbc_cli$breast_carcinoma_progesterone_receptor_status == 'Negative' &TCGA_tnbc_cli$lab_proc_her2_neu_immunohistochemistry_receptor_status == 'Negative', ]
TCGA_tnbc_cli$A0_Samples <- paste0(TCGA_tnbc_cli$A0_Samples, '-01')
TCGA_tnbc_cli <- TCGA_tnbc_cli[, c("A0_Samples", "A1_OS", "A2_Event")]
colnames(TCGA_tnbc_cli) <- c("A0_Samples", "OS.time", "OS")
table(TCGA_tnbc_cli$OS)
TCGA_tnbc_cli$OS <- ifelse(TCGA_tnbc_cli$OS == 'Alive', 0, 1)
rownames(TCGA_tnbc_cli) <- TCGA_tnbc_cli$A0_Samples
TCGA_tnbc_samples <- intersect(TCGA_tnbc_cli$A0_Samples, colnames(TCGA_tnbc))
TCGA_tnbc_cli <- TCGA_tnbc_cli[TCGA_tnbc_samples, ]





intersect(rownames(TCGA_tnbc), dHRGs)
setdiff(dHRGs, rownames(TCGA_tnbc))
intersect('CCN2', rownames(TCGA_tnbc))
intersect('ADGRG1', rownames(TCGA_tnbc))
intersect('TLNRD1', rownames(TCGA_tnbc))

TCGA_tnbc_dHRGs <- intersect(dHRGs, rownames(TCGA_tnbc))
intersect(rownames(TCGA_tnbc), TCGA_tnbc_dHRGs)




TCGA_tnbc_tpm <- TCGA_tnbc[, TCGA_tnbc_cli$A0_Samples]

dim(TCGA_tnbc_tpm)
TCGA_tnbc_tpm_log2 <- log2(TCGA_tnbc_tpm + 1)
boxplot(TCGA_tnbc_tpm_log2)






library(ggplot2)
TCGA_forest_log2 <- t(TCGA_tnbc_tpm_log2[, TCGA_tnbc_cli$A0_Samples])
TCGA_forest_log2 <- cbind(TCGA_tnbc_cli, TCGA_forest_log2)
TCGA_genes_index <- which(colnames(TCGA_forest_log2)%in%TCGA_tnbc_dHRGs)

TCGA_exp_symbol_processed.coxph <- cox_univariant_gene_regr(TCGA_forest_log2,TCGA_genes_index)
TCGA.forestplot.coxph <- TCGA_exp_symbol_processed.coxph
TCGA.forestplot.coxph <- TCGA.forestplot.coxph[TCGA.forestplot.coxph$Pvalue<0.05,]
dim(TCGA.forestplot.coxph)
TCGA.forestplot <- draw_forestplot(TCGA.forestplot.coxph)
TCGA.forestplot

TCGA.fp.coxph <- TCGA_exp_symbol_processed.coxph
TCGA.fp.table <- TCGA.fp.coxph[,c("gName","Pvalue","HR","Low.95.CI.","High.95.CI.")]
TCGA.fp.table <- TCGA.fp.table[order(TCGA.fp.table$HR,decreasing=F),]
genes_order <- as.character(TCGA.fp.table$gName)
rownames(TCGA.fp.coxph) <- TCGA.fp.coxph[,1]
genes_order_CI <- TCGA.fp.coxph[genes_order[length(genes_order):1],2:5]
genes_order_CI_sig <- genes_order_CI[which(genes_order_CI[,1]<0.05),]
write.table(genes_order_CI,"files/TCGA_genes_order_CI.txt",quote=F,row.names=T,col.names=T,sep="\t")
write.table(genes_order_CI_sig,"files/TCGA_genes_order_CI_p0.05.txt",quote=F,row.names=T,col.names=T,sep="\t")


intersect(as.character(TCGA.forestplot.coxph$gName), 
          as.character(GSE103091.forestplot.coxph$gName))
intersect(as.character(cbioportal.forestplot.coxph$gName), 
          as.character(TCGA.forestplot.coxph$gName))
intersect(as.character(cbioportal.forestplot.coxph$gName), 
          as.character(GSE103091.forestplot.coxph$gName))




comm_sig_genes <- intersect(as.character(TCGA.forestplot.coxph$gName), 
                            as.character(GSE103091.forestplot.coxph$gName))
comm_sig_genes

union_sig_genes <- unique(c(as.character(TCGA.forestplot.coxph$gName), 
                            as.character(GSE103091.forestplot.coxph$gName)))
union_sig_genes




forestplot <- ggpubr::ggarrange(GSE103091.forestplot,
                                TCGA.forestplot,
                                ncol = 2, nrow = 1,
                                labels = toupper(letters)[1:2],
                                align = "hv")   
forestplot
ggsave(plot = forestplot,
       filename = 'PDFs/univariateforestplot .pdf',
       width = 10, height = 6, device = cairo_pdf)


retrive_cluster_names <- function(myd,myd_consensusmap,hvalue){
  sample_names<-rownames(myd)
  myd_cut_list<-lapply(cut(myd_consensusmap$Colv,hvalue)$lower, function(l)rapply(l,function(i)i))
  cluster_sample_names<-c()
  tmp_cluster<-c()
  c_index<- 1
  for(i in myd_cut_list){
    cluster_sample_names<-c(cluster_sample_names,as.character(sample_names[unlist(i)]))
    tmp_cluster<-c(tmp_cluster,rep(paste("C",c_index,sep=""),length(unlist(i))))
    c_index<-c_index+1
  }
  cluster_df<-data.frame("Sample"=cluster_sample_names,"Cluster"=tmp_cluster)
  return(cluster_df);
}

setdiff(union_sig_genes, colnames(GSE103091_sig_cox_data))
GSE103091_union_sig_genes <- intersect(union_sig_genes, colnames(GSE103091_sig_cox_data))
library(NMF)
GSE103091_nmf_data <- t(GSE103091_sig_cox_data)
GSE103091_nmf <- nmf(GSE103091_nmf_data[GSE103091_union_sig_genes,], 2:10, nrun=50, seed=12345)
pdf('PDFs/GSE103091_NMF-1.pdf',width = 8,height = 6)
plot(GSE103091_nmf)
dev.off()

consensusmap(GSE103091_nmf,labCol=NA,labRow=NA,tracks=NA)



GSE103091_nmf_2 <- nmf(GSE103091_nmf_data[GSE103091_union_sig_genes,] , 2, nrun=50, seed=12345)
pdf('PDFs/GSE103091_NMF-2.pdf',width = 6,height = 6)
GSE103091_nmf_2_map <- consensusmap(GSE103091_nmf_2,
                                    labCol=NA,
                                    labRow=NA,
                                    tracks=NA)
dev.off()

GSE103091_nmf_2_cluster <- retrive_cluster_names(t(GSE103091_nmf_data),
                                                 GSE103091_nmf_2_map,
                                                 0.9)
table(GSE103091_nmf_2_cluster$Cluster)
colnames(GSE103091_nmf_2_cluster)<-c("Sample","Cluster")
rownames(GSE103091_nmf_2_cluster) <- GSE103091_nmf_2_cluster[,1]



GSE103091_cluster <- merge(GSE103091_cli, GSE103091_nmf_2_cluster, 
                          by.x = 'Samples', by.y = 'Sample')

GSE103091_NMF_OS <- data.frame(time = GSE103091_cluster$OS.time, 
                               OS = GSE103091_cluster$OS, 
                               Groups = GSE103091_cluster$Cluster)

GSE103091_NMF_km <- ggplotKMCox(GSE103091_NMF_OS,
                                labs = c('C1', 'C2'))
GSE103091_NMF_km


dat <- GSE103091_nmf_data[GSE103091_union_sig_genes, GSE103091_nmf_2_cluster$Sample]
annotation_col  <- data.frame(Groups = factor(GSE103091_nmf_2_cluster$Cluster))
rownames(annotation_col)=colnames(dat)

table(GSE103091_nmf_2_cluster$Cluster)
pdf('PDFs/GSE103091_nmf_heatmap.pdf',width = 6,height = 6)
bk=unique(c(seq(-2, 2, length=100)))
pheatmap(dat, scale = 'row', 
         breaks = bk,
         annotation_col = annotation_col,
         cluster_cols = F, cluster_rows = T,
         show_rownames = F, show_colnames = F,
         gaps_col = 34,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()


dat <- t(GSE103091_nmf_data[comm_sig_genes, GSE103091_nmf_2_cluster$Sample])
dat <- as.data.frame(dat)
groups <- GSE103091_nmf_2_cluster$Cluster
GSE103091_Violin <- groupViolin(dat, groups, ylab = 'log2(exp)', 
                          legend.pos = 'bl', xangle = 0)
GSE103091_Violin


library(NMF)
cbioportal_nmf_data <- cbioportal_exp[, cbioportal_tnbc_samples]

setdiff(union_sig_genes, rownames(cbioportal_nmf_data))
intersect('LMPHM4', rownames(cbioportal_nmf_data))
cbioportal_union_sig_genes <- c(intersect(union_sig_genes, rownames(cbioportal_nmf_data)))
cbioportal_nmf <- nmf(cbioportal_nmf_data[cbioportal_union_sig_genes,], 2:10, nrun=50, seed=12345)
pdf('PDFs/cbioportal_NMF-1.pdf',width = 8,height = 6)
plot(cbioportal_nmf)
dev.off()

consensusmap(cbioportal_nmf,labCol=NA,labRow=NA,tracks=NA)



cbioportal_nmf_2 <- nmf(cbioportal_nmf_data[cbioportal_union_sig_genes,], 
                        2, 
                        nrun=50, 
                        seed=12345)
pdf('PDFs/cbioportal_NMF-2.pdf',width = 6,height = 6)
cbioportal_nmf_2_map <- consensusmap(cbioportal_nmf_2,
                                    labCol=NA,
                                    labRow=NA,
                                    tracks=NA)
dev.off()

cbioportal_nmf_2_cluster <- retrive_cluster_names(t(cbioportal_nmf_data),
                                                 cbioportal_nmf_2_map,
                                                 0.9)
table(cbioportal_nmf_2_cluster$Cluster)


colnames(cbioportal_nmf_2_cluster)<-c("Sample","Cluster")
rownames(cbioportal_nmf_2_cluster) <- cbioportal_nmf_2_cluster[,1]
cbioportal_nmf_2_cluster <- arrange(cbioportal_nmf_2_cluster, Cluster)


cbioportal_exp_filter$samples <- rownames(cbioportal_exp_filter)
cbioportal_cluster <- merge(cbioportal_exp_filter, cbioportal_nmf_2_cluster, 
                           by.x = 'samples', by.y = 'Sample')

cbioportal_NMF_OS <- data.frame(time = cbioportal_cluster$OS.time, 
                               OS = cbioportal_cluster$OS, 
                               Groups = cbioportal_cluster$Cluster)

cbioportal_NMF_km <- ggplotKMCox(cbioportal_NMF_OS,
                                 labs = c('C1', 'C2'))
cbioportal_NMF_km


dat <- cbioportal_nmf_data[cbioportal_union_sig_genes, cbioportal_nmf_2_cluster$Sample]
annotation_col  <- data.frame(Groups = factor(cbioportal_nmf_2_cluster$Cluster))
rownames(annotation_col)=colnames(dat)


pdf('PDFs/cbioportal_nmf_heatmap.pdf',width = 6,height = 6)
bk=unique(c(seq(-1.5, 1.5, length=100)))
pheatmap(dat, scale = 'row', 
         breaks = bk,
         annotation_col = annotation_col,
         cluster_cols = F, cluster_rows = T,
         show_rownames = F, show_colnames = F,
         gaps_col = 138,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()


dat <- t(cbioportal_nmf_data[comm_sig_genes, cbioportal_nmf_2_cluster$Sample])
dat <- as.data.frame(dat)
groups <- cbioportal_nmf_2_cluster$Cluster
cbioportal_Violin <- groupViolin(dat, groups, ylab = 'log2(exp)', 
                                legend.pos = 'bl', xangle = 0)
cbioportal_Violin


library(NMF)
tcga_nmf_data <- TCGA_tnbc_tpm_log2

setdiff(union_sig_genes, rownames(tcga_nmf_data))
intersect('ADGRG1', rownames(tcga_nmf_data))
tcga_union_sig_genes <- c(intersect(union_sig_genes, rownames(tcga_nmf_data)))

tcga_nmf <- nmf(tcga_nmf_data[tcga_union_sig_genes,], 2:10, nrun=50, seed=12345)
pdf('PDFs/tcga_NMF-1.pdf',width = 8,height = 6)
plot(tcga_nmf)
dev.off()

consensusmap(tcga_nmf,labCol=NA,labRow=NA,tracks=NA)



tcga_nmf_2 <- nmf(tcga_nmf_data[tcga_union_sig_genes,], 
                        2, 
                        nrun=50, 
                        seed=12345)
pdf('PDFs/tcga_NMF-2.pdf',width = 6,height = 6)
tcga_nmf_2_map <- consensusmap(tcga_nmf_2,
                                     labCol=NA,
                                     labRow=NA,
                                     tracks=NA)
dev.off()

tcga_nmf_2_cluster <- retrive_cluster_names(t(tcga_nmf_data),
                                                  tcga_nmf_2_map,
                                                  0.9)
table(tcga_nmf_2_cluster$Cluster)
colnames(tcga_nmf_2_cluster)<-c("Sample","Cluster")
rownames(tcga_nmf_2_cluster) <- tcga_nmf_2_cluster[,1]


tcga_cluster <- merge(TCGA_tnbc_cli, tcga_nmf_2_cluster, 
                            by.x = 'A0_Samples', by.y = 'Sample')

tcga_NMF_OS <- data.frame(time = tcga_cluster$OS.time, 
                                OS = tcga_cluster$OS, 
                                Groups = tcga_cluster$Cluster)

tcga_NMF_km <- ggplotKMCox(tcga_NMF_OS,
                                 labs = c('C1', 'C2'))
tcga_NMF_km


dat <- tcga_nmf_data[tcga_union_sig_genes, tcga_nmf_2_cluster$Sample]
annotation_col  <- data.frame(Groups = factor(tcga_nmf_2_cluster$Cluster))
rownames(annotation_col)=colnames(dat)


pdf('PDFs/tcga_nmf_heatmap.pdf',width = 6,height = 6)
bk=unique(c(seq(-2, 2, length=100)))
pheatmap(dat, scale = 'row', 
         breaks = bk,
         annotation_col = annotation_col,
         cluster_cols = F, cluster_rows = T,
         show_rownames = F, show_colnames = F,
         gaps_col = 41,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()


dat <- t(tcga_nmf_data[comm_sig_genes, tcga_nmf_2_cluster$Sample])
dat <- as.data.frame(dat)
groups <- tcga_nmf_2_cluster$Cluster
tcga_Violin <- groupViolin(dat, groups, ylab = 'log2(exp)', 
                                 legend.pos = 'bl', xangle = 0)
tcga_Violin

NMF_km <- ggpubr::ggarrange(GSE103091_NMF_km,
                            tcga_NMF_km,
                            cbioportal_NMF_km,
                            ncol = 3,nrow = 1,
                            labels = toupper(letters)[1:3],
                            align = "hv")
NMF_km
ggsave(plot = NMF_km,
       filename = 'PDFs/NMF_km.pdf',
       width = 15, height = 5, device = cairo_pdf)

library(ggplot2)
data_Violin <- ggpubr::ggarrange(GSE103091_NMF_km,
                                 tcga_NMF_km,
                                 cbioportal_NMF_km,
                                 GSE103091_Violin,
                                 tcga_Violin,
                                 cbioportal_Violin,
                                 ncol = 3, nrow = 2,
                                 labels = toupper(letters)[4:9],
                                 align = "hv")
data_Violin
ggsave(plot = data_Violin,
       filename = 'PDFs/data_Violin.pdf',
       width = 15, height = 9, device = cairo_pdf)

save.image('TNBC-001.RData')


library(clusterRepro)
calculate_cluster_centroids <- function(expd,expd_cluster,colCentroid){
  unique(expd_cluster$Cluster)->clusters;
  if(colCentroid){
    t(expd)->expd;
  }
  c()->centroid_values;
  for(cl in clusters){
    which(expd_cluster$Cluster==cl)->cl_index;
    as.character(expd_cluster$Sample[cl_index])->cl_samples;
    expd[cl_samples,]->expd_cl;
    apply(expd_cl,2,function(x){mean(x)})->means;
    c(centroid_values,means)->centroid_values;
  }
  matrix(centroid_values,nrow=length(clusters),byrow=T)->centroid_values.matrix;
  clusters->rownames(centroid_values.matrix);
  colnames(expd)->colnames(centroid_values.matrix);
  return(centroid_values.matrix)
}


GSE103091_tmp <- GSE103091_sig_cox_data[, GSE103091_union_sig_genes]
GSE103091_tmp_cluster <- GSE103091_nmf_2_cluster
GSE103091_tmp_cluster$Cluster <- as.factor(GSE103091_tmp_cluster$Cluster)
GSE103091.centroids <- calculate_cluster_centroids(GSE103091_tmp,
                                                   GSE103091_tmp_cluster,
                                                   FALSE)
GSE103091.centroids <- t(GSE103091.centroids)
cbioportal_tmp <- cbioportal_nmf_data[cbioportal_union_sig_genes, cbioportal_nmf_2_cluster$Sample]
cbioportal_tmp <- as.matrix(cbioportal_tmp)
cbioportal_result1 <- clusterRepro(cbioportal_tmp, 
                                   GSE103091.centroids,
                                   Number.of.permutations = 1000)
cbioportal_result.IGP <- IGP.clusterRepro(GSE103091.centroids,
                                       Data=cbioportal_tmp)



tcga_tmp <- as.matrix(TCGA_tnbc_tpm_log2[tcga_union_sig_genes, ])
tcga_result1 <- clusterRepro(tcga_tmp, 
                                   GSE103091.centroids,
                                   Number.of.permutations = 1000)
tcga_result.IGP <- IGP.clusterRepro(GSE103091.centroids,
                                       Data=tcga_tmp)



library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(lumi)
tcga_methy <- fread('Z:/TCGA/Matrix/Methy/BRCA_450k', sep = '\t', header = T, 
                    stringsAsFactors = F, data.table = F)
rownames(tcga_methy) <- tcga_methy$sample
tcga_methy <- tcga_methy[, -1]
tcga_methy[1:5, 1:5]

tcga_methy_samples <- intersect(colnames(tcga_methy), colnames(tcga_tmp))
table(tcga_nmf_2_cluster[tcga_methy_samples, ]$Cluster)
tcga_methy <- tcga_methy[, tcga_methy_samples]
densityPlot(as.matrix(tcga_methy), main="Raw",
            legend=FALSE)

tcga_methy_nomr1 <- normalizeBetweenArrays(tcga_methy)

densityPlot(as.matrix(tcga_methy_nomr1), main="Normalized",
            legend=FALSE)

ann450kSub <- read.delim('files/methy_anno_filtered.txt', header = T,
                         stringsAsFactors = F)
table(ann450kSub$UCSC_RefGene_Group)


ann450kSub <- ann450kSub[ann450kSub[, 6] == 'Body' | ann450kSub[, 6] == 'TSS1500' | ann450kSub[, 6] == 'TSS200', ]

dHRGs

intersect(dHRGs, ann450kSub$UCSC_RefGene_Name)
methy_dHRGs <- intersect(dHRGs, ann450kSub$UCSC_RefGene_Name)


methy_dHRGs_index <- c()
for (gene in methy_dHRGs) {
  print(gene)
  ss <- which(ann450kSub$UCSC_RefGene_Name == gene)
  print(ss)
  methy_dHRGs_index <- c(methy_dHRGs_index, ss)
}
methy_dHRGs_re <- ann450kSub[methy_dHRGs_index, ]
rownames(methy_dHRGs_re) <- methy_dHRGs_re$Name

tcga_methy_cor_datas <- tcga_methy_nomr1[methy_dHRGs_re$Name, ]
tcga_methy_cor_datas <- t(na.omit(tcga_methy_cor_datas))
cortable <- methy_dHRGs_re[colnames(tcga_methy_cor_datas), ]

tcga_exp_cor_datas <- t(TCGA_tnbc_tpm_log2)
colnames(tcga_exp_cor_datas) <- gsub('CCN2', 'CTGF', colnames(tcga_exp_cor_datas))
colnames(tcga_exp_cor_datas) <- gsub('ADGRG1', 'GPR56', colnames(tcga_exp_cor_datas))
colnames(tcga_exp_cor_datas) <- gsub('TLNRD1', 'MESDC1', colnames(tcga_exp_cor_datas))
tcga_exp_cor_datas <- tcga_exp_cor_datas[rownames(tcga_methy_cor_datas),
                                         unique(cortable$UCSC_RefGene_Name)]

setdiff(cortable$UCSC_RefGene_Name, colnames(tcga_exp_cor_datas))


head(cortable)
tcga_methy_cor_datas[1:5, 1:5]
tcga_exp_cor_datas[1:5, 1:5]
dim(tcga_methy_cor_datas)
dim(tcga_exp_cor_datas)
dim(cortable)

cal_correlation_MET_EXP <- function(cortable, methydatas, expdatas){
  MetProbe <- c() 
  ExpGene <- c()  
  Cor_value <- c()      
  p_value <- c()   
  Groups <- c()
  for (i in 1:nrow(cortable)) {
    methy_name <- cortable[i, ]$Name
    gene_name <- cortable[i, ]$UCSC_RefGene_Name
    group <- cortable[i, ]$UCSC_RefGene_Group
    tmp_cor <- cor.test(methydatas[, methy_name], 
                        expdatas[, gene_name], 
                        method = 'spearman')
    MetProbe <- c(MetProbe, methy_name)
    ExpGene <- c(ExpGene, gene_name)
    Cor_value <- c(Cor_value, tmp_cor$estimate)
    p_value <- c(p_value, tmp_cor$p.value)
    Groups <- c(Groups, group)
  }
  res <- data.frame("MetProbes"   = MetProbe,
                    "ExpGenes"    = ExpGene,
                    "Cor"        = as.numeric(Cor_value),
                    "Cor.pvalue" = as.numeric(p_value),
                    "Groups"       = Groups)
  return(res)
}

correlation_MET_EXP <- cal_correlation_MET_EXP(cortable, 
                                               tcga_methy_cor_datas, 
                                               tcga_exp_cor_datas)
correlation_MET_EXP <- crbind2DataFrame(correlation_MET_EXP)
cor_MET_EXP_filter <- correlation_MET_EXP[correlation_MET_EXP$Cor.pvalue < 0.01 & abs(correlation_MET_EXP$Cor) > 0.6, ]
correlation_MET_EXP$groups <- ifelse(correlation_MET_EXP$Cor.pvalue < 0.01 & abs(correlation_MET_EXP$Cor) > 0.6, 'yes', 'no')
write.table(correlation_MET_EXP, file = 'files/correlation_MET_EXP.txt', sep = '\t',
            row.names = F, quote = F)

unique(cor_MET_EXP_filter$ExpGenes)


library(dplyr)
cor_MET_EXP_filter <- arrange(cor_MET_EXP_filter, ExpGenes)
tcga_methy_group <- tcga_nmf_2_cluster[rownames(tcga_methy_cor_datas), ]
tcga_methy_group <- arrange(tcga_methy_group, Cluster)
table(tcga_methy_group$Cluster)
labels_row <- cor_MET_EXP_filter$ExpGenes

annotation_col  <- data.frame(Groups = factor(tcga_methy_group$Cluster))
rownames(annotation_col) <- colnames(dat)

dat <- t(tcga_methy_cor_datas[tcga_methy_group$Sample, cor_MET_EXP_filter$MetProbes])

pdf('PDFs/TCGAmeth_gene.pdf', width = 6, height = 6)
bk=unique(c(seq(-2.5, 2.5, length=100)))
pheatmap::pheatmap(dat,
                   breaks = bk,
                   scale = 'row',
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                   cluster_rows = F,
                   cluster_cols = F,
                   gaps_col = 31,
                   show_colnames = F,
                   labels_row = labels_row,
                   border_color = NA,
                   annotation_col = annotation_col)
dev.off()

labels_row <- cor_MET_EXP_filter$MetProbes
pdf('PDFs/TCGAmeth_meth.pdf', width = 6, height = 6)
bk=unique(c(seq(-2.5, 2.5, length=100)))
pheatmap::pheatmap(dat,
                   breaks = bk,
                   scale = 'row',
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                   cluster_rows = F,
                   cluster_cols = F,
                   gaps_col = 31,
                   show_colnames = F,
                   labels_row = labels_row,
                   border_color = NA,
                   annotation_col = annotation_col)
dev.off()


correlation_MET_EXP$Cor.pvalue <- correlation_MET_EXP$Cor.pvalue + 1e-12
correlation_MET_EXP$log10 <- -log10(correlation_MET_EXP$Cor.pvalue)

methy_Quadrant_plot <- ggplot(correlation_MET_EXP, aes(x = Cor, y = log10, colour = groups)) +
  geom_point() + xlab("Correlation") + ylab("-log10(P.value)") + xlim(-1, 1) + ylim(0, 18) +
  geom_hline(aes(yintercept=-log10(0.01)), colour="#990000", linetype="dashed")+
  geom_vline(aes(xintercept=0.6), colour="#990000", linetype="dashed")+ 
  geom_vline(aes(xintercept=-0.6), colour="#990000", linetype="dashed") +
  scale_colour_manual(values=c(yes = "red", no = "grey"))
methy_Quadrant_plot

showText_methy <- correlation_MET_EXP[abs(correlation_MET_EXP$Cor) > 0.6 & correlation_MET_EXP$Cor.pvalue < 0.01, ]$ExpGenes
methy_Quadrant_plot <- mg_volcano(correlation_MET_EXP$Cor, 
           correlation_MET_EXP$Cor.pvalue,
           cutFC=0.6,
           cutPvalue=0.01,
           legend.pos='tr',
           ylab='-log10(P)',
           xlab='Correlation')
methy_Quadrant_plot
ggsave(plot = methy_Quadrant_plot,
       filename = 'PDFs/methy_Quadrant_plot.pdf',
       width = 5, height = 5, device = cairo_pdf)

library(corrplot)
methy_cor_genes <- unique(cor_MET_EXP_filter$ExpGenes)

m6a <- read.delim('files/m6a-genes.txt', header = T, stringsAsFactors = F)
m6a_genes <- m6a$Symbols

tcga_m6a_datas <- TCGA_tnbc_tpm_log2
rownames(tcga_m6a_datas) <- gsub('RBMX', 'HNRNPG', rownames(tcga_m6a_datas))
rownames(tcga_m6a_datas) <- gsub('VIRMA', 'KIAA1429', rownames(tcga_m6a_datas))



tcga_m6a_datas <- t(tcga_m6a_datas[c(methy_cor_genes, m6a_genes), ])
tcga_m6a_cor <- cor(tcga_m6a_datas, method = 'spearman')
tcga_m6a_cor_ci <- cor.mtest(tcga_m6a_datas,
                          conf.level=0.95)

change_cor_p2symbols<-function(cor_testd,cutps){
  cor_testd$p->cor_pd;
  c()->cor_p_symbols;
  for(i in 1:nrow(cor_pd)){
    for(j in 1:ncol(cor_pd)){
      if(cor_pd[i,j]<cutps[1]){
        c(cor_p_symbols,"***")->cor_p_symbols;
      }else if(cor_pd[i,j]<cutps[2] && cor_pd[i,j]>=cutps[1]){
        c(cor_p_symbols,"**")->cor_p_symbols;
      }else if(cor_pd[i,j]<cutps[3] && cor_pd[i,j]>=cutps[2]){
        c(cor_p_symbols,"*")->cor_p_symbols;
      }else{
        c(cor_p_symbols," ")->cor_p_symbols;
      }
    }
  }
  matrix(cor_p_symbols,ncol=ncol(cor_pd),byrow=T)->cor_p_symbols;
  return(cor_p_symbols);
}
m6a_genes_p2symbols <- change_cor_p2symbols(tcga_m6a_cor_ci, c(0.001,0.01,0.05))

rownames(m6a_genes_p2symbols) <- colnames(tcga_m6a_cor)
colnames(m6a_genes_p2symbols) <- colnames(tcga_m6a_cor)


prepare_text_colors<-function(catd,genes){
  unique(catd$m6a_process)->catd_types;
  mg_colors[1:length(catd_types)]->catd_types_col;
  catd_types->names(catd_types_col);
  #--
  as.character(catd$m6a_process)->m6a_process;
  catd$Symbol->names(m6a_process);
  m6a_process[genes]->m6a_process;
  catd_types_col[m6a_process]->catd_colors
  return(catd_colors)
}
m6a_process_types_colors <- prepare_text_colors(m6a ,m6a_genes)

library(superheat)
pdf('PDFs/m6a_dHRGS_cor.pdf', width = 8, height = 6)
superheat(tcga_m6a_cor[methy_cor_genes, m6a_genes],
          X.text=m6a_genes_p2symbols[methy_cor_genes, m6a_genes],
          bottom.label.text.angle = 90,
          bottom.label.text.size = 4,
          heat.pal = colorRampPalette(c("navy", "white", "firebrick3"))(100),
          bottom.label.text.col=m6a_process_types_colors)
dev.off()
write.table(tcga_m6a_cor[methy_cor_genes, m6a_genes], 
            file = 'files/tcga_m6a_cor.txt', sep = '\t',
            quote = F)
write.table(m6a_genes_p2symbols[methy_cor_genes, m6a_genes], 
            file = 'files/m6a_genes_p2symbols.txt', sep = '\t',
            quote = F)



library(data.table)
load('origin_datas/TCGA/CNV/brca_cnv_seg_data.RData')
rownames(brca_cnv_seg_data) <- brca_cnv_seg_data$Sample
com_dHRGs_cnv_genes <- intersect(colnames(brca_cnv_seg_data), dHRGs)
TCBC_cnv_seg_data <- brca_cnv_seg_data[, com_dHRGs_cnv_genes]
TCBC_cnv_seg_data <- as.data.frame(TCBC_cnv_seg_data)
TCBC_cnv_cn_data <- (2**TCBC_cnv_seg_data) * 2

tcga_exp_cnv_cor_datas <- t(TCGA_tnbc_tpm_log2)
tcga_exp_cnv_cor_datas <- tcga_exp_cnv_cor_datas[tcga_cluster$A0_Samples, com_dHRGs_cnv_genes]


cal_correlation_CNV_EXP <- function(CNVdatas, expdatas){
  ExpGene <- c()  
  Cor_value <- c()    
  p_value <- c()  
  genes <- intersect(colnames(CNVdatas), colnames(expdatas))
  samples <- intersect(rownames(CNVdatas), rownames(expdatas))
  CNVdatas <- CNVdatas[samples, genes]
  expdatas <- expdatas[samples, genes]
  
  for (gene in genes) {
    tmp_cor <- cor.test(CNVdatas[, gene], 
                        expdatas[, gene], 
                        method = 'spearman')
    ExpGene <- c(ExpGene, gene)
    Cor_value <- c(Cor_value, tmp_cor$estimate)
    p_value <- c(p_value, tmp_cor$p.value)
  }
  res <- data.frame("ExpGenes"    = ExpGene,
                    "Cor"        = as.numeric(Cor_value),
                    "Cor.pvalue" = as.numeric(p_value))
  return(res)
}

correlation_CNV_EXP <- cal_correlation_CNV_EXP(TCBC_cnv_cn_data, 
                                               tcga_exp_cnv_cor_datas)
correlation_CNV_EXP <- crbind2DataFrame(correlation_CNV_EXP)
cor_CNV_EXP_filter <- correlation_CNV_EXP[correlation_CNV_EXP$Cor.pvalue < 0.01 & abs(correlation_CNV_EXP$Cor) > 0.5, ]
correlation_CNV_EXP$groups <- ifelse(correlation_CNV_EXP$Cor.pvalue < 0.01 & abs(correlation_CNV_EXP$Cor) > 0.5, 'yes', 'no')

write.table(correlation_CNV_EXP, 
            file = 'files/correlation_CNV_EXP.txt', sep = '\t',
            quote = F, row.names = F)

unique(cor_CNV_EXP_filter$ExpGenes)


library(dplyr)
cor_CNV_EXP_filter <- arrange(cor_CNV_EXP_filter, ExpGenes)


tcga_cnv_group <- tcga_nmf_2_cluster[rownames(TCBC_cnv_cn_data), ]
tcga_cnv_group <- arrange(tcga_cnv_group, Cluster)
table(tcga_cnv_group$Cluster)
labels_row <- cor_CNV_EXP_filter$ExpGenes
annotation_col  <- data.frame(Groups = factor(tcga_cnv_group$Cluster))
rownames(annotation_col) <- tcga_cnv_group$Sample

dat <- t(TCBC_cnv_cn_data[tcga_cnv_group$Sample, cor_CNV_EXP_filter$ExpGenes])

pdf('PDFs/TCGACNV_gene.pdf', width = 6, height = 6)
bk=unique(c(seq(-2.5, 2.5, length=100)))
pheatmap::pheatmap(dat,
                   breaks = bk,
                   scale = 'row',
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                   cluster_rows = T,
                   cluster_cols = F,
                   gaps_col = 41,
                   show_colnames = F,
                   labels_row = labels_row,
                   border_color = NA,
                   annotation_col = annotation_col)
dev.off()


correlation_CNV_EXP$Cor.pvalue <- correlation_CNV_EXP$Cor.pvalue + 1e-12


showText_cnv <- correlation_CNV_EXP[abs(correlation_CNV_EXP$Cor) > 0.5 & correlation_CNV_EXP$Cor.pvalue < 0.01, ]$ExpGenes
CNV_Quadrant_plot <- mg_volcano(correlation_CNV_EXP$Cor, 
                                  correlation_CNV_EXP$Cor.pvalue,
                                  cutFC=0.5,
                                  cutPvalue=0.01,
                                  legend.pos='tl',
                                  ylab='-log10(P)',
                                  xlab='Correlation')
CNV_Quadrant_plot
ggsave(plot = CNV_Quadrant_plot,
       filename = 'PDFs/CNV_Quadrant_plot.pdf',
       width = 6, height = 4, device = cairo_pdf)




library(maftools)
brca_mutect <- read.maf('Z://TCGA/Matrix/mutect2/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.gz')
oncoplot(maf = brca_mutect, top = 10, bgCol = "#FFFFFF")
brca_mutect_cli <- brca_mutect@clinical.data
brca_mutect_cli <- crbind2DataFrame(brca_mutect_cli)
class(brca_mutect_cli)
brca_mutect_cli$samples <- substr(brca_mutect_cli$Tumor_Sample_Barcode, 1,15)

brca_mutect_cli <- merge(brca_mutect_cli, tcga_nmf_2_cluster, 
                         by.x = 'samples', by.y = "Sample")
brca_mutect_cli$samples <- substr(brca_mutect_cli$samples, 1, 12)

brac_cli <- read.delim('origin_datas/TCGA/Clinical BCR XML.merge.txt', header = T,
                       stringsAsFactors = F)
brac_cli <- merge(brac_cli, brca_mutect_cli, by.x = 'A0_Samples', by.y = 'samples')
write.table(brac_cli, file = 'origin_datas/TCGA/Clinical BCR XML.merge.filter.txt',
            sep = '\t', row.names = F, quote = F)


brca_mutect <- read.maf(maf = 'Z://TCGA/Matrix/mutect2/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.gz', clinicalData = 'origin_datas/TCGA/Clinical BCR XML.merge.filter.txt', verbose = FALSE)
pdf('PDFs/TCGA-BRCA-mutect-plot.pdf', width = 12, height = 8)
oncoplot(maf = brca_mutect,
         # colors = vc_cols,
         bgCol = "#FFFFFF",
         clinicalFeatures = c('Cluster', 'A2_Event'),
         sortByAnnotation = T,
         genes = c('TP53', 'PIK3CA', 'RB1', 'PTEN', 'NF1', 'BRAF','BRCA1')
         # top = 10
         )
dev.off()


maf.cli.samples <- substr(getClinicalData(brca_mutect)$Tumor_Sample_Barcode, 1, 15)
com_samples <- intersect(maf.cli.samples, tcga_cluster$A0_Samples)

PIK3CA <- as.character(subsetMaf(maf = brca_mutect, genes = c('PIK3CA'), mafObj = FALSE)$Tumor_Sample_Barcode)
PIK3CA <- substr(PIK3CA, 1, 15)
PIK3CA_samples <- intersect(PIK3CA, com_samples)
NO_PIK3CA_samples <- setdiff(com_samples, PIK3CA_samples)
PIK3CA_Mutant <- c(rep('YES', length(PIK3CA_samples)), rep('NO', length(NO_PIK3CA_samples)))
PIK3CA_Mutant <- data.frame(cbind(sample = c(PIK3CA_samples, NO_PIK3CA_samples), PIK3CA_Mutant))
tcga_cluster_mut <- merge(tcga_cluster, PIK3CA_Mutant, 
                      by.x = 'A0_Samples', by.y = 'sample', all = TRUE)


TP53 <- as.character(subsetMaf(maf = brca_mutect, genes = c('TP53'), mafObj = FALSE)$Tumor_Sample_Barcode)
TP53 <- substr(TP53, 1, 15)
TP53_samples <- intersect(TP53, com_samples)
NO_TP53_samples <- setdiff(com_samples, TP53_samples)
TP53_Mutant <- c(rep('YES', length(TP53_samples)), rep('NO', length(NO_TP53_samples)))
TP53_Mutant <- data.frame(cbind(sample = c(TP53_samples, NO_TP53_samples), TP53_Mutant))
tcga_cluster_mut <- merge(tcga_cluster_mut, TP53_Mutant, 
                          by.x = 'A0_Samples', by.y = 'sample', all = TRUE)

PTEN <- as.character(subsetMaf(maf = brca_mutect, genes = c('PTEN'), mafObj = FALSE)$Tumor_Sample_Barcode)
PTEN <- substr(PTEN, 1, 15)
PTEN_samples <- intersect(PTEN, com_samples)
NO_PTEN_samples <- setdiff(com_samples, PTEN_samples)
PTEN_Mutant <- c(rep('YES', length(PTEN_samples)), rep('NO', length(NO_PTEN_samples)))
PTEN_Mutant <- data.frame(cbind(sample = c(PTEN_samples, NO_PTEN_samples), PTEN_Mutant))
tcga_cluster_mut <- merge(tcga_cluster_mut, PTEN_Mutant, 
                          by.x = 'A0_Samples', by.y = 'sample', all = TRUE)

NF1 <- as.character(subsetMaf(maf = brca_mutect, genes = c('NF1'), mafObj = FALSE)$Tumor_Sample_Barcode)
NF1 <- substr(NF1, 1, 15)
NF1_samples <- intersect(NF1, com_samples)
NO_NF1_samples <- setdiff(com_samples, NF1_samples)
NF1_Mutant <- c(rep('YES', length(NF1_samples)), rep('NO', length(NO_NF1_samples)))
NF1_Mutant <- data.frame(cbind(sample = c(NF1_samples, NO_NF1_samples), NF1_Mutant))
tcga_cluster_mut <- merge(tcga_cluster_mut, NF1_Mutant, 
                          by.x = 'A0_Samples', by.y = 'sample', all = TRUE)

RB1 <- as.character(subsetMaf(maf = brca_mutect, genes = c('RB1'), mafObj = FALSE)$Tumor_Sample_Barcode)
RB1 <- substr(RB1, 1, 15)
RB1_samples <- intersect(RB1, com_samples)
NO_RB1_samples <- setdiff(com_samples, RB1_samples)
RB1_Mutant <- c(rep('YES', length(RB1_samples)), rep('NO', length(NO_RB1_samples)))
RB1_Mutant <- data.frame(cbind(sample = c(RB1_samples, NO_RB1_samples), RB1_Mutant))
tcga_cluster_mut <- merge(tcga_cluster_mut, RB1_Mutant, 
                          by.x = 'A0_Samples', by.y = 'sample', all = TRUE)

EGFR <- as.character(subsetMaf(maf = brca_mutect, genes = c('EGFR'), mafObj = FALSE)$Tumor_Sample_Barcode)
EGFR <- substr(EGFR, 1, 15)
EGFR_samples <- intersect(EGFR, com_samples)
NO_EGFR_samples <- setdiff(com_samples, EGFR_samples)
EGFR_Mutant <- c(rep('YES', length(EGFR_samples)), rep('NO', length(NO_EGFR_samples)))
EGFR_Mutant <- data.frame(cbind(sample = c(EGFR_samples, NO_EGFR_samples), EGFR_Mutant))
tcga_cluster_mut <- merge(tcga_cluster_mut, EGFR_Mutant, 
                          by.x = 'A0_Samples', by.y = 'sample', all = TRUE)

BRAF <- as.character(subsetMaf(maf = brca_mutect, genes = c('BRAF'), mafObj = FALSE)$Tumor_Sample_Barcode)
BRAF <- substr(BRAF, 1, 15)
BRAF_samples <- intersect(BRAF, com_samples)
NO_BRAF_samples <- setdiff(com_samples, BRAF_samples)
BRAF_Mutant <- c(rep('YES', length(BRAF_samples)), rep('NO', length(NO_BRAF_samples)))
BRAF_Mutant <- data.frame(cbind(sample = c(BRAF_samples, NO_BRAF_samples), BRAF_Mutant))
tcga_cluster_mut <- merge(tcga_cluster_mut, BRAF_Mutant, 
                          by.x = 'A0_Samples', by.y = 'sample', all = TRUE)

table(tcga_cluster_mut[tcga_cluster_mut$Cluster == 'C1', ]$BRAF_Mutant)
table(tcga_cluster_mut[tcga_cluster_mut$Cluster == 'C2', ]$BRAF_Mutant)

tdata <- matrix(c(0,36,2,62), ncol = 2)
tdata
chisq.test(tdata)
fisher.test(tdata)





cbioportal_roc_km_datas <- cbioportal_cluster
tcga_roc_km_datas <- as.data.frame(t(TCGA_tnbc_tpm_log2[TCGA_tnbc_dHRGs, ]))
tcga_roc_km_datas$samples <- rownames(tcga_roc_km_datas)
tcga_roc_km_datas <- merge(tcga_cluster, tcga_roc_km_datas,
                           by.x = "A0_Samples", by.y = "samples")

GSE103091_roc_km_datas <- GSE103091_exp_log2[, GSE103091_dHRGs]
GSE103091_roc_km_datas <- crbind2DataFrame(GSE103091_roc_km_datas)
GSE103091_roc_km_datas$samples <- rownames(GSE103091_roc_km_datas)
GSE103091_roc_km_datas <- merge(GSE103091_cli, GSE103091_roc_km_datas, 
                                by.x = "Samples", by.y = "samples")

roc_km_datas <- GSE103091_roc_km_datas

library(survival)
fmla <- as.formula(paste0("Surv(OS.time, OS) ~"
                          ,paste0(comm_sig_genes,collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(roc_km_datas))
lan <- coef(cox)
lan
genes <- names(cox$coefficients)
genes


cbioportal.tr <- as.numeric(lan%*%as.matrix(t(cbioportal_roc_km_datas[,genes])))

cbioportal.tr <- mosaic::zscore(cbioportal.tr)
cutoff <- 0

cbioportal.km <- ggplotKMCox(data.frame(cbioportal_roc_km_datas$OS.time,
                       cbioportal_roc_km_datas$OS,
                       ifelse(cbioportal.tr >= cutoff,'H','L')),
                       title = 'METABRIC',
                       labs = c('High', 'Low'))
cbioportal.km
cbioportal.roc <- ggplotTimeROC(cbioportal_roc_km_datas$OS.time,
              cbioportal_roc_km_datas$OS,
              cbioportal.tr,
              mks = c(1,3,5))
cbioportal.roc

cbioportal.km.roc <- plotCoxModel_Batch(cbioportal.tr,
                   cbioportal_roc_km_datas[,genes],
                   cbioportal_roc_km_datas$OS.time,
                   cbioportal_roc_km_datas$OS,
                   cutoff = 0)
cbioportal.km.roc
ggsave(plot = cbioportal.km.roc,
       filename = 'cbioportal.km.roc.pdf',
       width = 13, height = 10, device = cairo_pdf)

tcga.tr <- as.numeric(lan%*%as.matrix(t(tcga_roc_km_datas[,genes])))

tcga.tr <- mosaic::zscore(tcga.tr)
cutoff <- 0


tcga.tr <- tcga.tr - 0.15

tcga.km <- ggplotKMCox(data.frame(tcga_roc_km_datas$OS.time,
                       tcga_roc_km_datas$OS,
                       ifelse(tcga.tr >= 0,'H','L')),
                       title = 'TCGA',
            labs = c('High', 'Low'))
tcga.km
tcga.roc <- ggplotTimeROC(tcga_roc_km_datas$OS.time,
              tcga_roc_km_datas$OS,
              tcga.tr,
              mks = c(1,3,5))
tcga.roc

tcga.km.roc <- plotCoxModel_Batch(tcga.tr,
                                        tcga_roc_km_datas[,genes],
                                        tcga_roc_km_datas$OS.time,
                                        tcga_roc_km_datas$OS,
                                        cutoff = 0)
tcga.km.roc
ggsave(plot = tcga.km.roc,
       filename = 'tcga.km.roc.pdf',
       width = 13, height = 10, device = cairo_pdf)

GSE103091.tr <- as.numeric(lan%*%as.matrix(t(GSE103091_roc_km_datas[,genes])))
cutoff <- median(GSE103091.tr)
GSE103091.tr <- mosaic::zscore(GSE103091.tr)
cutoff <- 0


fit <- survivalROC::survivalROC(Stime = GSE103091_roc_km_datas$OS.time, 
                   status = GSE103091_roc_km_datas$OS, 
                   marker = GSE103091.tr, 
                   predict.time = 365*1, 
                   method = "KM")
fit$AUC
cutoff <- fit$cut.values[which.max(fit$TP - fit$FP)]
cutoff
GSE103091.tr <- GSE103091.tr - cutoff

GSE103091.km <- ggplotKMCox(data.frame(GSE103091_roc_km_datas$OS.time,
                       GSE103091_roc_km_datas$OS,
                       ifelse(GSE103091.tr >= 0,'H','L')),
                       title = 'GSE103091',
            labs = c('High', 'Low'))
GSE103091.km
GSE103091.roc <- ggplotTimeROC(GSE103091_roc_km_datas$OS.time,
              GSE103091_roc_km_datas$OS,
              GSE103091.tr,
              mks = c(2,3,5))
GSE103091.roc

GSE103091.km.roc <- plotCoxModel_Batch(GSE103091.tr,
                                  GSE103091_roc_km_datas[,genes],
                                  GSE103091_roc_km_datas$OS.time,
                                  GSE103091_roc_km_datas$OS,
                                  cutoff = 0,
                                  mks = c(2,3,5))
GSE103091.km.roc
ggsave(plot = GSE103091.km.roc,
       filename = 'GSE103091.km.roc.pdf',
       width = 13, height = 10, device = cairo_pdf)


datas_km_roc <- ggpubr::ggarrange(tcga.roc,
                                  GSE103091.roc,
                                  cbioportal.roc,
                                  tcga.km,
                                  GSE103091.km,
                                  cbioportal.km,
                                  labels = toupper(letters)[1:6],
                                  align = "hv")

datas_km_roc
ggsave(plot = datas_km_roc,
       filename = 'PDFs/datas_km_roc.pdf',
       width = 15, height = 10, device = cairo_pdf)

ggsave(plot = datas_km_roc,
       filename = 'datas_km_roc.pdf',
       width = 15, height = 10, device = cairo_pdf)

