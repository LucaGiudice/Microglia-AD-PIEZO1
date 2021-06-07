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
op1_data_path="output/op1_data_rdy.rda"
markers_path="data/mmu_markers.txt"
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
markers=fread(markers_path,data.table = F)
#multicoreParam <- MulticoreParam(workers = 16)

### Step 1 ----
if(S1 & S2){
  sce <- as.SingleCellExperiment(seuObj)
  sce = cxds_bcds_hybrid(sce, verb=TRUE)
  sce_woDup=sce[,sce$hybrid_score<0.35]
  seuObj <- as.Seurat(sce_woDup)
  
  #Compute percentage of mt
  seuObj[["percent.mt"]] <- PercentageFeatureSet(seuObj, pattern = "^MT-")
  
  #Quality control and filtering of cells based on counts, features and mt
  VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  #Filtering
  seuObj <- subset(seuObj, subset = nCount_RNA > 300 & nCount_RNA < 9000 & 
                     nFeature_RNA > 300 & nFeature_RNA < 5600)
  
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
  n_pca=30
  seuObj <- ScoreJackStraw(seuObj, dims = 1:n_pca)
  JackStrawPlot(seuObj, dims = 1:n_pca)
  #Good quality of the data, I use all the PCA components
  n_best_dims=30
  seuObj <- FindNeighbors(seuObj, dims = 1:n_best_dims)
  #I decided to minimize small difficult to annotate clusters
  seuObj <- FindClusters(seuObj, resolution = 0.06)
  head(Idents(seuObj), 5)
  seuObj <- RunUMAP(seuObj, dims = 1:n_best_dims)
  
  #png("./output/op1_umap.png", units="in", width=22, height=12, res=300)
    DimPlot(seuObj, reduction = "umap", label=T, pt.size = 0.8, label.size = 8)
  #dev.off()
  
    #Filter to keep the biggest 10 clusters
    seuObj_filt=seuObj[,as.numeric((as.character(Idents(seuObj))))<11]
    #Clustering and UMAP plot
    n_best_dims=30
    seuObj_filt <- FindNeighbors(seuObj_filt, dims = 1:n_best_dims)
    seuObj_filt <- FindClusters(seuObj_filt, resolution = 0.06)
    seuObj_filt <- RunUMAP(seuObj_filt, dims = 1:n_best_dims)
    DimPlot(seuObj_filt, reduction = "umap", label=T, pt.size = 0.8, label.size = 8)
    
  #Preserve this final big step to replicate in the future the same analysis
  save(seuObj,seuObj_filt,file=step2_lock)
  
}else{
  load(step2_lock)
}

types=unique(markers$Cell_type)
top_markers_l=list()
for(type in types){
  gs=markers$Genes[markers$Cell_type==type]
  top_markers_l[[type]]=gs
}

#Set data
geno=as.matrix(seuObj_filt@assays[["RNA"]]@data[,])
original_i=Idents(object = seuObj_filt)

#Perform GSVA marker test
geno2=t(geno)
geno2=aggregate(geno2, list(original_i), mean)
names=geno2[,1];geno2=geno2[,-1];
geno2=t(geno2);colnames(geno2)=names

#Perform GSVA score analysis https://f1000research.com/articles/8-296
gsva_res=gsva(geno2,gset.idx.list=top_markers_l,method="gsva",abs.ranking=F,
              min.sz=1,max.sz=999,parallel.sz=16,mx.diff=T);gsva_res=t(gsva_res)
ann_v="start"
for(k_cl in 1:nrow(gsva_res)){
  ann_v=c(ann_v,colnames(gsva_res)[which.max(gsva_res[k_cl,])])
}
ann_v=ann_v[-1]
gsva_ann=data.frame(cls=rownames(gsva_res),ann=ann_v,stringsAsFactors = F)

#Update the annotation of the clusters based on scCATCH result
orig_clusters=Idents(object = seuObj_filt)
lvs=levels(orig_clusters);names(lvs)=lvs
lvs[gsva_ann[,1]]=as.character(gsva_ann$ann)
seuObj_filt <- RenameIdents(seuObj_filt, lvs)

svg(filename=out_plot_path1, width=15, height=15, pointsize=12)
DimPlot(seuObj_filt, reduction = "umap", label=T, 
        pt.size = 0.8, label.size = 8)
dev.off()

svg(filename=out_plot_path2, width=15, height=15, pointsize=12)
FeaturePlot(seuObj_filt, features = c("Piezo1"), pt.size = 2.5,
            min.cutoff = 1, max.cutoff = 3, cols=c("grey85","darkorange1"),
            order=T, label=T, label.size = 6)
dev.off()

png(out_plot_path3, units="in", width=22, height=12, res=300)
  p1=VlnPlot(seuObj_filt, features = c("Piezo1"))
  p2=VlnPlot(seuObj_filt, features = c("Piezo1"), split.by = "group")
  p1+p2
dev.off()

condXP = FetchData(object = seuObj_filt, 
                   vars = c('ident','treatment',"model","group",'Piezo1'))
save(condXP,file="./output/op2_PIEZO1_condXP.rda")

info=seuObj_filt@meta.data
info$cls=Idents(object = seuObj_filt)
save(geno,info,seuObj_filt,file="output/op2_data_ready.rda")

seuObj_filt@meta.data$orig.ident=Idents(object = seuObj_filt)
seuObj_filt@meta.data$ident=Idents(object = seuObj_filt)






