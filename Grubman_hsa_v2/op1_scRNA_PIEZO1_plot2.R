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
op_scRNA_analysis="./output/op1_subidXP.rda"
#Output
out_boxplot="output/op1/op1_submic_PIEZO1_exprs.png"
out_excel="output/op1/op1_submic_PIEZO1_exprs.xlsx"

#Load and normalize data ----
cat("Load data \n")
load(op_scRNA_analysis)

#Set data
x=subidXP

#Find the names of the microglia subgroups X condition
#Remove the ZZZ groups
subgs=unique(x$subIDmXcond)[c(-1,-6)]
subgs=subgs[order(subgs)]

#Filter out rows not belonging to the subgroups
keep_indxs=0
orig_subgs=x$subIDmXcond
for(subg in subgs){
  keep_indxs=c(keep_indxs,grep(subg,orig_subgs,fixed = T))
}
keep_indxs=keep_indxs[-1]
x=x[keep_indxs,]
x=x[,-2]

#Divide the groups
subIDm=substr(x$subIDmXcond,1,2)
condition=substr(x$subIDmXcond,4,6)
x=cbind(x,subIDm,condition)
colnames(x)[c(2,3,4)]=c("Expression_value","Microglia_subgroup","Condition")
x$Condition=gsub("ct","CONTROL",x$Condition)

#Plot
x=x[x$Expression_value!=0,]
p2=ggplot(x, aes(x=Microglia_subgroup, y=Expression_value, fill=Condition)) + geom_boxplot()
p2=p2 + theme(text = element_text(size=30)) + scale_fill_manual(values=c("gold1", "grey"))
p2
write.xlsx(x, out_excel, sheetName="boxplot", append=F, row.names = F)

png(out_boxplot, units="in", width=22, height=12, res=300)
  plot(p2)
dev.off()



