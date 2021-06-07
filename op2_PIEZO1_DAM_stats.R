#Clean workspace and memory ----
#https://support.bioconductor.org/p/84439/
rm(list=ls())
gc()

#Set working directory----
gps0=getwd()
gps0=paste(gps0,"/%s",sep="")
rootDir=gps0
setwd(gsub("%s","",rootDir))
set.seed(8)

#Load libraries----

library("data.table")
library("matrixStats")
library("plyr")
library("purrr")
library("factoextra")
library("statmod")
library("matrixStats")
library("ggplot2")
library("ggpubr")
library("xlsx")
library("Seurat")
library("ggrepel")
see = function(m,nr=5,nc=5){
  if(nrow(m)<5){
    nr=nrow(m)
  }
  if(ncol(m)<5){
    nc=ncol(m)
  }
  
  m[1:nr,1:nc]
}
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
get_res_DAM = function(matrix_ra_freq,matrix_exp,thr=0.7){
  matrix_freq=cbind(matrix_ra_freq,rowMins(matrix_ra_freq[,1:3]))
  colnames(matrix_freq)[ncol(matrix_freq)]="min_freq"
  matrix_freq=matrix_freq[matrix_freq[,ncol(matrix_freq)]>=thr,]
  matrix_exp=matrix_exp[rownames(matrix_freq),]
  m=matrix_exp[abs(rowSums(matrix_exp))!=1,]
  cor_m=cor(t(m[,1:3]))
  m_top=m[abs(cor_m[,1])>=thr,]
  top_cor_m=cor(t(m_top[,1:3]))
  top_cor_df=data.frame(Genes=rownames(top_cor_m),Cor=top_cor_m[,1],stringsAsFactors = F)
  matrix_freq_df=as.data.frame(matrix_freq)
  matrix_freq_df$Genes=rownames(matrix_freq_df)
  res_DAM=merge(top_cor_df,matrix_freq_df,by="Genes")
  return(res_DAM)
}

#Set variables ----
#Input
op1_stats2="./Henna_scRNA_Shaul_mmu_v2/output/op3_DAM_analysis.rda"
op1_stats3="./Henna_scRNA_Zhou_mmu/output/op5_DAM_analysis.rda"
#Output
out_boxplot1="output_summary/op2_mmu_PIEZO1.png"
out_excel1="output_summary/op2_mmu_PIEZO1_freqs.xlsx"
out_excel2="output_summary/op2_mmu_PIEZO1_expre.xlsx"
out_rda="output_summary/op2_mmu_PIEZO1_stats.rda"

#Load and normalize data ----
cat("Load data \n")
dss_names=c("Shaul_cor_DAM","Zhou_cor_DAM")
load(op1_stats2)
Shaul_DAM=get_res_DAM(matrix_ra_freq,matrix_exp,thr = 0.7)
Shaul_DAM_all=get_res_DAM(matrix_ra_freq,matrix_exp,thr = 0.1)
load(op1_stats3)
Zhou_DAM=get_res_DAM(matrix_ra_freq,matrix_exp,thr = 0.7)
Zhou_DAM_all=get_res_DAM(matrix_ra_freq,matrix_exp,thr = 0.1)

DAM_profiles=merge(Shaul_DAM,Zhou_DAM,by="Genes")
DAM_profiles=DAM_profiles[(DAM_profiles$Cor.x>0 & DAM_profiles$Cor.y>0)|(DAM_profiles$Cor.x<0 & DAM_profiles$Cor.y<0),]
DAM_profiles_all=merge(Shaul_DAM_all,Zhou_DAM_all,by="Genes")
DAM_profiles_all=DAM_profiles_all[(DAM_profiles_all$Cor.x>0 & DAM_profiles_all$Cor.y>0)|(DAM_profiles_all$Cor.x<0 & DAM_profiles_all$Cor.y<0),]
DAM_profiles_all$names=""
DAM_profiles_all$names[match(DAM_profiles$Genes,DAM_profiles_all$Genes)]=DAM_profiles$Genes

#http://www.sthda.com/english/articles/32-r-graphics-essentials/131-plot-two-continuous-variables-scatter-graph-and-alternatives/
#https://www.mzes.uni-mannheim.de/socialsciencedatalab/article/datavis/
df=DAM_profiles_all
df$min_cells=rowMins(as.matrix(df[,c(7,13)]))*100
df$DAMs_labels.x <- as.factor(df$DAMs_labels.x)
.labs <- df$names
svg(filename="./output_summary/PIEZO1_DAM_corr_analysis.svg", width=15, height=15, pointsize=12)
ggplot(df, aes(x = Cor.x, y = Cor.y, size = min_cells)) + geom_point(shape = 18) +
  geom_point(aes(color = DAMs_labels.x)) +
  geom_label_repel(aes(label = .labs,  color = DAMs_labels.x), size = 5) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme(text = element_text(size=20)) +
  xlab("Correlation DAM-PIEZO1 in Shaul et al.") + 
  ylab("Correlation DAM-PIEZO1 in Zhou et al.") +
  labs(size = "Frequency in PIEZO1 cells") +
  labs(color = "DAM down - PIEZO1 - DAM up")
dev.off()


colnames(df)=gsub("\\.x",".Shaul",colnames(df))
colnames(df)=gsub("\\.y",".Zhou",colnames(df))
write.xlsx(df,file="./output_summary/PIEZO1_DAM_corr_analysis.xlsx",sheetName = "PIEZO1_DAM_cor_mmu",
           col.names = T, row.names = F, append=F)
up_cor=df[df$Cor.Shaul>0 & df$Cor.Zhou>0 & nchar(df$names)>1,]
dw_cor=df[df$Cor.Shaul<0 & df$Cor.Zhou<0 & nchar(df$names)>1,]
write.xlsx(up_cor,file="./output_summary/PIEZO1_DAM_corr_analysis.xlsx",sheetName = "up_DAM_cor",
           col.names = T, row.names = F, append=T)
write.xlsx(dw_cor,file="./output_summary/PIEZO1_DAM_corr_analysis.xlsx",sheetName = "dw_DAM_cor",
           col.names = T, row.names = F, append=T)












