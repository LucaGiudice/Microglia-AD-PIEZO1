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
library("xlsx")
library("Rfast")

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
gsm_path="data/GSE98969_RAW/"
gsm_files=list.files(path=gsm_path,all.files = TRUE,full.names = TRUE)
gsm_files=gsm_files[-seq(1,2)]
gsm_names=list.files(path=gsm_path,all.files = TRUE,full.names = FALSE)
gsm_names=gsm_names[-seq(1,2)]
output_op1="output/op1_data_rdy.rda"

#Load metadata ---
meta_df = fread("data/GSE98969_experimental_design_f.txt",data.table = F)
#Remove useless info
keep=c(seq(1,3),5,7)
meta_df = meta_df[,keep]
#Add a column related to condition AD or WT
meta_df$Cond=substr(meta_df$Batch_desc,1,2)
#Build df with only AD and WT
vec=substr(meta_df$Batch_desc,1,4)
AD_i=grep(pattern="AD\\dm", x=vec)
WT_i=grep(pattern="WT\\dm", x=vec)
meta_df_AD = meta_df[AD_i,]
del=-grep(unique(meta_df_AD$Batch_desc)[7],meta_df_AD$Batch_desc,fixed = T)
meta_df_AD = meta_df_AD[del,]
meta_df_WD = meta_df[WT_i,]
#Finish it
meta_df = rbind(meta_df_AD,meta_df_WD)

#Load and assemble files ----
for(k in 1:length(gsm_files)){
  cat(k,"on",length(gsm_files),"\n")
  t1=fread(gsm_files[k],data.table = FALSE)
  rownames(t1)=t1[,1];t1=t1[,-1];
  if(k==1){
    scData=t1
  }else{
    scData=cbind(scData,t1)
  }
}

#Processing ----
#remove zero rows
keep_bol=(rowsums(as.matrix(scData),parallel = TRUE))!=0
scData=scData[keep_bol,]

#Remove cells not having meta info
keep_cells=!is.na(match(colnames(scData),meta_df$Well_ID))
scData=scData[,keep_cells]
keep_cells=!is.na(match(meta_df$Well_ID,colnames(scData)))
meta_df=meta_df[keep_cells,]
keep_cells=match(meta_df$Well_ID,colnames(scData))
keep_cells=keep_cells[!is.na(keep_cells)]
scData=scData[,keep_cells]
rownames(meta_df)=meta_df$Well_ID

#Output
save(scData,meta_df,file=output_op1)



