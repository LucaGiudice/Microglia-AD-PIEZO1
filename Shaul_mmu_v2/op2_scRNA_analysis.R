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
library("GSVA")

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
op1_data_path="output/op1_data_rdy.rda"
markers_path="data/mmc1.xlsx"
de_mic_path="data/mmc2.xlsx"
#Step locks
step1_lock="keep/op2_seuObj_s1.rda";S1=F;
step2_lock="keep/op2_seuObj_s2.rda";S2=F;if(S2){S1=T};
step2_ann="keep/op2_clu_annotation.rda";
#output
output_dir="output/op2"
out_plot_path1="./output/op2/op2_umap_ann.svg"
out_plot_path2="./output/op2/op2_umap_PIEZO1.svg"
out_plot_path3="./output/op2/op2_exp_PIEZO1.svg"
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}else{
  print("Dir already exists!")
}

#Load data ----
load(op1_data_path)
markers=read_xlsx(markers_path)
de_mic=read_excel(de_mic_path)

#Check sample order in geno and info ---
ord=seq(1,nrow(meta_df))
matches=match(colnames(scData),meta_df$Well_ID)
cat("sample order in geno and info is identical:",all.equal(matches,ord),"\n")

### Step 1 ----
if(S1 & S2){
  seuObj <- CreateSeuratObject(counts = scData, 
                             project = "HK", 
                             min.cells = 20, 
                             min.features = 200,
                             meta.data = meta_df)
  
  #Compute percentage of mt
  seuObj[["percent.mt"]] <- PercentageFeatureSet(seuObj, pattern = "^MT-")
  
  #Quality control and filtering of cells based on counts, features and mt
  #VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  #Filtering
  seuObj <- subset(seuObj, subset = nCount_RNA < 5000 & nFeature_RNA < 2000)
  
  #Normalization
  seuObj <- NormalizeData(seuObj)
  
  #Feature Selection
  seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)
  
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(seuObj), 10)
  
  #Scaling
  all.genes <- rownames(seuObj)
  seuObj <- ScaleData(seuObj, features = all.genes)
  
  #PCA
  seuObj <- RunPCA(seuObj, features = VariableFeatures(object = seuObj))
  #Clustering
  seuObj <- JackStraw(seuObj, num.replicate = 100, dims=30)
  #Save
  save(seuObj,file=step1_lock)
}

### Step 2 ----
if(S2){
  if(!S1){
    load(step1_lock)
  }
  #JackStraw to determine best projections on low dimensional space
  n_pca=30
  seuObj <- ScoreJackStraw(seuObj, dims = 1:n_pca)
  JackStrawPlot(seuObj, dims = 1:n_pca)
  #Good quality of the data, I use all the PCA components
  n_best_dims=20
  seuObj <- FindNeighbors(seuObj, dims = 1:n_best_dims)
  #I decided to minimize small difficult to annotate clusters
  seuObj <- FindClusters(seuObj, resolution = 0.2)
  head(Idents(seuObj), 5)
  seuObj <- RunUMAP(seuObj, dims = 1:n_best_dims)

  #png("./output/op1_umap.png", units="in", width=22, height=12, res=300)
    DimPlot(seuObj, reduction = "umap", label=T, 
            pt.size = 0.8, label.size = 8)
  #dev.off()
  #Preserve this final big step to replicate in the future the same analysis
  save(seuObj,file=step2_lock)
  
}else{
  load(step2_lock)
  load(step2_ann)
}

#Get the top expressed of each cell type in the study
markers=markers[,-2]
markers=as.data.frame(markers)
rownames(markers)=markers[,1];markers=markers[,-1]
markers=as.matrix(markers)

#Generate df and list
genes=rownames(markers)
top_markers=apply(markers,2,function(x){
  y=genes[order(x,decreasing = T)]
  z=y[1:10]
  return(z)
})
top_markers_l=list()
for(i in 1:ncol(top_markers)) {
  top_markers_l[[i]] <- top_markers[ , i]
}
names(top_markers_l)=colnames(top_markers)

#Set microglia 2 and 3 markers from the DE table
mark_mic2=c("Apoe","Csf1r","Sparc","Cx3cr1","Serinc3","P2ry12","Malat1")
de_mic3vs1=de_mic[,c(1,2)]
colnames(de_mic3vs1)=c("genes","lFC")
de_mic3vs1=as.data.frame(de_mic3vs1)
mark3=de_mic3vs1$genes[order(de_mic3vs1$lFC,decreasing = T)]
top_mark3=mark3[1:15]
top_markers_l[["Microglia2."]]=mark_mic2
top_markers_l[["Microglia3."]]=top_mark3

#Set data
geno=as.matrix(seuObj@assays[["RNA"]][,])
original_i=Idents(object = seuObj)

#Perform GSVA marker test
geno2=t(geno)
geno2=aggregate(geno2, list(original_i), mean)
names=geno2[,1];geno2=geno2[,-1];
geno2=t(geno2);colnames(geno2)=names

#Perform GSVA score analysis https://f1000research.com/articles/8-296
gsva_res=gsva(geno2,gset.idx.list=top_markers_l,method="gsva",abs.ranking=F,
              min.sz=5,max.sz=999,parallel.sz=16,mx.diff=T);gsva_res=t(gsva_res)
ann_v="start"
for(k_cl in 1:nrow(gsva_res)){
  ann_v=c(ann_v,colnames(gsva_res)[which.max(gsva_res[k_cl,])])
}
ann_v=ann_v[-1]
gsva_ann=data.frame(cls=rownames(gsva_res),ann=ann_v)

#Perform manual analysis on the genes
cls=as.character(unique(Idents(object = seuObj)))
cls=cls[order(cls,decreasing = F)]
cl_ann_v="start"
cl_ann_mic3_v=cl_ann_mic2_v=0
for(cl in cls){
  indxs=which(cl==original_i)
  geno_cl=geno[,indxs]
  gs_exp=rowMeans(geno_cl)
  
  ma_exp_v=0
  for(k_mark in 1:ncol(top_markers)){
    marks=top_markers[,k_mark]
    ma_exp=mean(gs_exp[marks])
    ma_exp_v=c(ma_exp_v,ma_exp)
  }

  cl_ann_mic3_v=c(cl_ann_mic3_v,mean(gs_exp[top_mark3]))
  cl_ann_mic2_v=c(cl_ann_mic2_v,mean(gs_exp[mark_mic2]))
  
  ma_exp_v=ma_exp_v[-1]
  cl_ann=colnames(top_markers)[which.max(ma_exp_v)]
  cl_ann_v=c(cl_ann_v,cl_ann)
}
cl_ann_v=cl_ann_v[-1]
df_ann=data.frame(orig_cls=cls,ann=cl_ann_v,mic2=cl_ann_mic2_v[-1],mic3=cl_ann_mic3_v[-1])
df_ann$ann=sub("\\,.*", "", df_ann$ann)
df_ann[4,2]="Microglia3"
df_ann[3,2]="Microglia2"
df_ann[7,2]="Mature B-cells"
df_ann[10,2]="Immature B-cells"

#Update the annotation of the clusters based on scCATCH result
orig_clusters=Idents(object = seuObj)
lvs=levels(orig_clusters);names(lvs)=lvs
lvs[df_ann[,1]]=df_ann$ann
names(lvs) <- levels(seuObj)
seuObj <- RenameIdents(seuObj, lvs)

svg(filename=out_plot_path1, width=15, height=15, pointsize=12)
DimPlot(seuObj, reduction = "umap", label=T, 
        pt.size = 0.8, label.size = 8)
dev.off()

svg(filename=out_plot_path2, width=15, height=15, pointsize=12)
FeaturePlot(seuObj, features = c("Fam38a"), pt.size = 2.5,
            min.cutoff = 1, max.cutoff = 3, cols=c("grey85","darkorange1"),
            order=T, label=T, label.size = 6)
dev.off()

png(out_plot_path3, units="in", width=22, height=12, res=300)
  p1=VlnPlot(seuObj, features = c("Fam38a"), split.by = "Cond", split.plot = TRUE,
             cols=c("yellow","grey"))
  p2=VlnPlot(seuObj, features = c("Fam38a"))
  p1+p2
dev.off()

condXP = FetchData(object = seuObj, vars = c('ident', 'Cond','Fam38a'))
save(condXP,file="./output/op2_PIEZO1_condXP.rda")

info=seuObj@meta.data
info$cls=Idents(object = seuObj)
save(geno,info,seuObj,file="output/op2_data_ready.rda")









