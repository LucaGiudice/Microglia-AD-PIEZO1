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

#Set variables ----
#Input
op1_stats1="./Henna_scRNA_Grubman_hsa_v2/output/op1/op1_hsa_PIEZO1_stats.rda"
op1_stats2="./Henna_scRNA_Shaul_mmu_v2/output/op2/op2_mmu_PIEZO1_stats.rda"
op1_stats3="./Henna_scRNA_Zhou_mmu/output/op2/op2_mmu_PIEZO1_stats.rda"
op1_stats4="./Henna_scRNA_Zhou_mmu/output/op4/op4_mmu_PIEZO1_stats.rda"
#Output
out_boxplot1="output_summary/op2_mmu_PIEZO1.png"
out_excel1="output_summary/op2_mmu_PIEZO1_freqs.xlsx"
out_excel2="output_summary/op2_mmu_PIEZO1_expre.xlsx"
out_rda="output_summary/op2_mmu_PIEZO1_stats.rda"

#Load and normalize data ----
cat("Load data \n")
dss_names=c("Grubman_hsa","Shaul_mmu","Zhou_mmu","Zhou_mmu_submicroglia")
load(op1_stats1)
stats1=cbind(info_df,exp_df)
load(op1_stats2)
stats2=cbind(info_df,exp_df)
load(op1_stats3)
stats3=cbind(info_df,exp_df)
load(op1_stats4)
stats4=cbind(info_df,exp_df)
dss_l=list(stats1,stats2,stats3,stats4)


app=F
for(k in 1:4){
  dss_name=dss_names[k]
  
  stats=dss_l[[k]]
  stats[,1]=range01(stats[,1])
  stats[,2]=range01(stats[,2])
  stats[,7]=range01(stats[,7])
  stats=stats[,-4]
  stats=stats[,c(3,4,7)]
  
  s_dss_name=paste("PIEZO1_stats",dss_name,sep="_")
  write.xlsx(stats,file="./output_summary/PIEZO1_stats.xlsx",sheetName = s_dss_name,
             col.names = T, row.names = F, append=app);app=T;
  
  stats_m=melt(stats)
  plot_s_dss_name=paste(s_dss_name,".png",sep="")
  p1 <- ggplot(stats_m,aes(cell_types,value,fill=variable)) + 
    geom_bar(stat="sum",na.rm=TRUE) + theme(text = element_text(size=20))
  png(plot_s_dss_name, units="in", width=22, height=12, res=300)
  plot(p1)
  dev.off()
}
