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
library("readxl")
library("GSVA")
library("DropletUtils")
library("BiocParallel")
library("scds")

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
op1_data_path="output/op2_data_ready.rda"
markers_path="data/mmu_markers.txt"
#Step locks
step1_lock="keep/op4_seuObj_s1.rda";S1=F;
step2_lock="keep/op4_seuObj_s2.rda";S2=F;if(S2){S1=T};
step2_ann="keep/op4_clu_annotation.rda";
#output
output_dir="output/op4"
out_plot_path1="./output/op4/op4_umap_ann.svg"
out_plot_path2="./output/op4/op4_umap_PIEZO1.svg"
out_plot_path3="./output/op4/op4_exp_PIEZO1.svg"
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}else{
  print("Dir already exists!")
}

#Load data ----
load(op1_data_path)
markers=fread(markers_path,data.table = F)
#multicoreParam <- MulticoreParam(workers = 16)
seuObj=subset(x = seuObj_filt, idents = "Microglia")
rm(seuObj_filt,geno,info)
### Step 1 ----
if(S1 & S2){
  #Compute percentage of mt
  seuObj[["percent.mt"]] <- PercentageFeatureSet(seuObj, pattern = "^MT-")
  
  #Quality control and filtering of cells based on counts, features and mt
  VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  #Filtering
  seuObj <- subset(seuObj, subset = nCount_RNA < 2000 & nFeature_RNA < 1500)
  
  #Normalization
  seuObj <- NormalizeData(seuObj)
  
  #Feature Selection
  seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 3000)
  
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
  n_pca=20
  seuObj <- ScoreJackStraw(seuObj, dims = 1:n_pca)
  JackStrawPlot(seuObj, dims = 1:n_pca)
  #Good quality of the data, I use all the PCA components
  n_best_dims=16
  seuObj <- FindNeighbors(seuObj, dims = 1:n_best_dims)
  #I decided to minimize small difficult to annotate clusters
  seuObj <- FindClusters(seuObj, resolution = 0.05)
  head(Idents(seuObj), 5)
  seuObj <- RunUMAP(seuObj, dims = 1:n_best_dims)

  ann_cls=c("0","WT_5XFAD","2","3")
  seuObj <- RenameIdents(object = seuObj, '0' = ann_cls[1], "1"=ann_cls[2], 
                          "2"=ann_cls[3], "3"=ann_cls[4])
  seuObj@meta.data$ann_cls=Idents(object = seuObj)
  #png("./output/op1_umap.png", units="in", width=22, height=12, res=300)
    DimPlot(seuObj, reduction = "umap", label=T, pt.size = 0.8, label.size = 8)
  #dev.off()
    
  png(out_plot_path3, units="in", width=22, height=12, res=300)
    VlnPlot(seuObj, features = c("Lpcat2","Tmem119","Ank","Piezo1","P2ry12","Csf1"))
  dev.off()
    
  #Preserve this final big step to replicate in the future the same analysis
  save(seuObj,file=step2_lock)
  
}else{
  load(step2_lock)
}

condXP = FetchData(object = seuObj, vars = c('ann_cls',"Lpcat2","Tmem119","Ank","Piezo1","P2ry12","Csf1"))
save(condXP,file="./output/op4_PIEZO1_condXP.rda")

info=seuObj@meta.data
info$cls=Idents(object = seuObj)
geno=as.matrix(seuObj@assays[["RNA"]]@data[,])
save(geno,info,seuObj,file="output/op4_data_ready.rda")



