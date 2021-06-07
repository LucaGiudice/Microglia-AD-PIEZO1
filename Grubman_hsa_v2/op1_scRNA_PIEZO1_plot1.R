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


#Set variables ----
#Input
op1_scRNA_analysis="./output/op1_data_ready.rda"
op1_scRNA_fetchData="./output/op1_batchesXP.rda"
#Output
out_boxplot1="output/op1/op1_hsa_PIEZO1.png"
out_excel1="output/op1/op1_hsa_PIEZO1_freqs.xlsx"
out_excel2="output/op1/op1_hsa_PIEZO1_expre.xlsx"
out_rda="output/op1/op1_hsa_PIEZO1_stats.rda"

#Load and normalize data ----
cat("Load data \n")
load(op1_scRNA_analysis)
load(op1_scRNA_fetchData)

#Set data
m=as.matrix(geno[,])
cell_types=as.character(info$cls)

#Analysis preview: Piezo1 in MUSC
classes=unique(cell_types)
classes=classes[-c(3,11)]
for(cl_k in 1:length(classes)){
  cl=classes[cl_k]
  cl_indxs=which(cell_types==cl)
  
  #Eval the percentage of expression of each gene
  cl_m=m[,cl_indxs]
  gs_perc_v=apply(cl_m,1,function(x){
    g_perc=(sum(x!=0)*100)/length(x)
    return(g_perc)
  })
  gs_n_v=apply(cl_m,1,function(x){
    g_n=sum(x!=0)
    return(g_n)
  })
  #Compute the rank of each gene perc
  gs_rank_v=rank(gs_perc_v,ties.method = "min")/length(gs_perc_v)
  #Get Piezo1 info
  Piezo1_n_cells=gs_n_v["PIEZO1"]
  Piezo1_perc=gs_perc_v["PIEZO1"]
  Piezo1_rank=gs_rank_v["PIEZO1"]
  Piezo1_df=data.frame(Piezo1_n_cells=Piezo1_n_cells,
                       Piezo1_perc=Piezo1_perc,
                       Piezo1_rank_perc=Piezo1_rank,
                       cell_types=cl,stringsAsFactors = F)
  rownames(Piezo1_df)=cl
  
  #Eval how much Piezo1 is expressed in its cells
  Piezo1_dist=cl_m["PIEZO1",]
  cells_indxs=which(Piezo1_dist!=0)
  if(length(cells_indxs)>2){
    Piezo1_m=cl_m[,cells_indxs]
    
    Piezo1_rank_m=apply(Piezo1_m,2,function(x){
      y=rank(x,ties.method ="min")
      z=y/length(y)
      return(z)
    })
    rownames(Piezo1_rank_m)=rownames(cl_m)
    Piezo1_rank_exp=Piezo1_rank_m["PIEZO1",]
    Piezo1_rank_exp=mean(Piezo1_rank_exp)
    Piezo1_exp=mean(Piezo1_dist[cells_indxs])
  }else{
    Piezo1_rank_exp=0
    Piezo1_exp=0
  }

  #Eval how much Piezo1 is expressed in all cells
  cl_m_NA=cl_m
  cl_m_NA[cl_m_NA==0]=NA
  type_exp=rowMeans(cl_m_NA,na.rm = T)
  Piezo1_rank_exp_all=rank(type_exp,ties.method ="min")
  Piezo1_rank_exp_all=Piezo1_rank_exp_all/length(Piezo1_rank_exp_all)
  Piezo1_rank_exp_all=Piezo1_rank_exp_all["PIEZO1"]

  Piezo1_rank_exp_df=data.frame(Piezo1_rank_exp_Pcells=Piezo1_rank_exp,
                                Piezo1_rank_exp_all=Piezo1_rank_exp_all,
                                Piezo1_exp_Pcells=Piezo1_exp,
                                cell_types=cl,
                                stringsAsFactors = F)
  rownames(Piezo1_rank_exp_df)=cl
  
  if(cl_k==1){
    info_df=Piezo1_df
    exp_df=Piezo1_rank_exp_df
  }else{
    info_df=rbind(info_df,Piezo1_df)
    exp_df=rbind(exp_df,Piezo1_rank_exp_df)
  }
}
save(info_df,exp_df,file=out_rda)

png("./output/op1/op1_perc_PIEZO1.png", units="in", width=22, height=12, res=300)
info_df=info_df[,c(1,3,4)]
ggplot(data=info_df, aes(x=cell_types, y=Piezo1_rank_perc)) + geom_bar(stat="identity") +
  theme(text = element_text(size=20)) +
  scale_y_continuous(breaks=seq(0,1,by=0.10))
dev.off()

write.xlsx(info_df,file="./output/op1/op1_perc_PIEZO1.xlsx",sheetName = "perc_PIEZO1",
           col.names = T, row.names = F)

png("./output/op1/op1_exp_PIEZO1.png", units="in", width=22, height=12, res=300)
exp_df=exp_df[,c(2,4)]
ggplot(data=exp_df, aes(x=cell_types, y=Piezo1_rank_exp_all)) + geom_bar(stat="identity") +
  theme(text = element_text(size=20)) +
  scale_y_continuous(breaks=seq(0,1,by=0.10))
dev.off()
write.xlsx(exp_df,file="./output/op1/op1_exp_PIEZO1.xlsx",sheetName = "exp_PIEZO1",
           col.names = T, row.names = F)

#Remove zero values
x2=batchesXP
rownames(x2)=seq(1,nrow(x2))
colnames(x2)=c("Condition","Cell_types","Expression_value")
x2$Condition=gsub("ct","CONTROL",x2$Condition)
x2$Condition=as.character(x2$Condition)
x2$Cell_types=as.character(x2$Cell_types)
x2=x2[x2$Expression_value!=0,]

#Plot how much Piezo1 is expressed in the clusters
p2=ggplot(x2, aes(x=Cell_types, y=Expression_value, fill=Condition)) + geom_boxplot()
p2=p2 + theme(text = element_text(size=18)) + scale_fill_manual(values=c("gold1", "grey"))

#Save the table of the plot
write.xlsx(x2, "./output/op1/op1_raw_exp_PIEZO1xcond.xlsx", sheetName="boxplot", 
           append=T, row.names = F)

png("./output/op1/op1_raw_exp_PIEZO1xcond.png", units="in", width=22, height=12, res=300)
  plot(p2)
dev.off()






