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
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140511
#https://www.nature.com/articles/s41591-019-0695-9#data-availability
gsm_path="data/GSE140511_RAW/"
gsm_files=list.files(path=gsm_path,all.files = TRUE,full.names = TRUE)
gsm_files=gsm_files[-seq(1,2)]
gsm_names=list.files(path=gsm_path,all.files = TRUE,full.names = FALSE)
gsm_names=gsm_names[-seq(1,2)]
output_op1="output/op1_data_rdy.rda"

#Create metadata ---
gsm_names=substr(gsm_names,12,50)
gsm_names=gsub("_barcodes.tsv","",gsm_names)
gsm_names=gsub("_features.tsv","",gsm_names)
gsm_names=gsub("_matrix.mtx","",gsm_names)
gsm_names=gsub("_genes.tsv","",gsm_names)
gsm_names_l=strsplit(gsm_names,"_",fixed = T)
meta_df=data.frame(matrix("CO",length(gsm_names),2),stringsAsFactors = F)
colnames(meta_df)=c("treatment","model")
for(k_l in 1:nrow(meta_df)){
  row=gsm_names_l[[k_l]]
  row=row[-length(row)]
  if(length(row)==1){
    meta_df[k_l,1]=row[1]
  }
  if(length(row)==2){
    meta_df[k_l,1:2]=row[1:2]
  }
  if(length(row)==3){
    row[2]=paste(row[2],row[3],sep="_")
    meta_df[k_l,1:2]=row[1:2]
  }
}

#Group name for samples -----
meta_df$group=paste(meta_df$treatment,meta_df$model,sep="_")


#Load and assemble files ----
mtx_indxs=grep("_matrix.mtx",gsm_files)
for(mtx_indx in mtx_indxs){
  cat(mtx_indx,"\n")
  mtx_file=gsm_files[mtx_indx]
  fea_file=gsm_files[mtx_indx-1]
  bar_file=gsm_files[mtx_indx-2]
  mtx=Matrix::readMM(mtx_file)
  genes=fread(fea_file,header = F,data.table = F)
  rownames(mtx)=genes$V2
  
  name_samples=meta_df$group[mtx_indx]
  colnames(mtx)=rep(name_samples,ncol(mtx))
  
  meta_df2=meta_df[rep(mtx_indx,ncol(mtx)),]
  if(mtx_indx==3){
    meta_dfs=meta_df2
    final_mtx=mtx
  }else{
    meta_dfs=rbind(meta_dfs,meta_df2)
    final_mtx=merge(final_mtx,mtx,by="row.names")
    rownames(final_mtx)=final_mtx[,1];final_mtx=final_mtx[,-1]
  }
}
colnames(final_mtx)=make.unique(meta_dfs$group)
rownames(meta_dfs)=colnames(final_mtx)
rm(mtx,meta_df2);gc()

seuObj <- CreateSeuratObject(counts = final_mtx, 
                             project = "HK", 
                             min.cells = 20, 
                             min.features = 200,
                             meta.data = meta_dfs)

#Output
save(seuObj,file=output_op1)



