#Clean workspace and memory ----
rm(list=ls())
gc()

#Set working directory----
gps0=getwd()
gps0=paste(gps0,"/%s",sep="")
rootDir=gps0
setwd(gsub("%s","",rootDir))
set.seed(8)

#Load libraries----
library('Seurat')
library("scater")
library("data.table")
library("matrixStats")
library("plyr")
library("doParallel")
library("parallel")
library("readxl")
library("Rfast")
library("scCATCH")
library("Rfast2")
library("Rfast")
library("ggplot2")
library("xlsx")

see = function(m,nr=5,nc=5){
  if(nrow(m)<5){
    nr=nrow(m)
  }
  if(ncol(m)<5){
    nc=ncol(m)
  }
  
  m[1:nr,1:nc]
}

#Set variables ----
#Input
DAM_markers_path="../Henna_scRNA_Shaul_mmu_v2/output/op3_DAM_analysis.rda"
geno_path="output/op4_data_ready.rda"
#Output
output_dir="output/op5"
out_DAM_markers="output/op5/op5_DAM_markersWPiezo1.xlsx"
out_boxplot1="output/op5/op5_DAM_markers_pos.png"
out_boxplot2="output/op5/op5_DAM_markers_freq.png"
out_DAM_analysis="output/op5_DAM_analysis.rda"
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}else{
  print("Dir already exists!")
}

#Load data ----
load(DAM_markers_path)
load(geno_path)

#Retrieve the DAM markers ----
DAM_up=rownames(matrix_exp)[matrix_exp[,4]==1]
DAM_dw=rownames(matrix_exp)[matrix_exp[,4]==-1]
DAMs=c(DAM_up,DAM_dw)
DAMs_labels=c(rep(1,length(DAM_up)),rep(-1,length(DAM_dw)))
keep=DAMs %in% rownames(geno)
DAMs=DAMs[keep]
DAMs_labels=c(0,DAMs_labels[keep])

#Cell types of interest
types=c("0","2","WT_5XFAD")
matrix_freq=matrix(0,length(DAMs)+1,3)
colnames(matrix_freq)=types;rownames(matrix_freq)=c("Piezo1",DAMs)
matrix_freq=cbind(matrix_freq,DAMs_labels)
matrix_ra_exp=matrix_ra_freq=matrix_exp=matrix_freq
for(k_type in 1:length(types)){

  #Retrieve the cluster cells in which PIEZO1 is expressed
  bidentx=which(info$cls==types[k_type])
  scData=geno[,bidentx]
  indx=match("Piezo1",rownames(scData))
  exb=scData[indx,]!=0
  Piezo1_cells=scData[,exb]
  #Analyse each gene for ratio of frequency and expression
      
      gs_perc_v=apply(Piezo1_cells,1,function(x){
        g_perc=((sum(x!=0))*100)/length(x)
        return(g_perc)
      })
      #Compute the rank of each gene perc
      gs_rank_perc_v=rank(gs_perc_v,ties.method = "min")/length(gs_perc_v)
      #Get genes frequency ratio
      Piezo1_ra_perc=gs_rank_perc_v["Piezo1"]
      DAM_ra_perc=gs_rank_perc_v[DAMs]
      matrix_ra_freq[,k_type]=c(Piezo1_ra_perc,DAM_ra_perc)
      
      #Get genes frequency
      Piezo1_perc=gs_perc_v["Piezo1"]
      DAM_perc=gs_perc_v[DAMs]
      matrix_freq[,k_type]=c(Piezo1_perc,DAM_perc)
    
      #Compute the rank of each gene expression
      gs_exp_v=apply(Piezo1_cells,1,function(x){
        x=x[x!=0]
        g_avg_exp=mean(x)
        if(is.na(g_avg_exp)){
          g_avg_exp=0
        }
        return(g_avg_exp)
      })
      #Compute the rank of each gene exp
      gs_rank_exp_v=rank(gs_exp_v,ties.method = "min")/length(gs_exp_v)
      #Get genes expression ratio info
      Piezo1_ra_exp=gs_rank_exp_v["Piezo1"]
      DAM_ra_exp=gs_rank_exp_v[DAMs]
      matrix_ra_exp[,k_type]=c(Piezo1_ra_exp,DAM_ra_exp)
      
      #Get genes expression
      Piezo1_exp=gs_exp_v["Piezo1"]
      DAM_exp=gs_exp_v[DAMs]
      matrix_exp[,k_type]=c(Piezo1_exp,DAM_exp)
}
matrix_ra_freq=round(matrix_ra_freq,5)
matrix_freq=round(matrix_freq,5)
matrix_ra_exp=round(matrix_ra_exp,5)
matrix_exp=round(matrix_exp,5)
# 
# m=matrix_exp
# m=m[abs(rowSums(m))!=1,]
# cor_m=cor(t(m[,1:3]))
# m_top=m[abs(cor_m[,1])>=0.9,]
# top_cor_m=cor(t(m_top[,1:3]))

save(matrix_ra_freq,matrix_freq,matrix_ra_exp,matrix_exp,file=out_DAM_analysis)

