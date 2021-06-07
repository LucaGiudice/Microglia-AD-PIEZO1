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
#https://www.nature.com/articles/s41593-019-0539-4#data-availability
#http://adsn.ddnetbio.com/
#http://adsn.ddnetbio.com/#tab-4691-1
#Input
geno_path="data/scRNA_rawCounts.tsv"
header_path="data/scRNA_metadata_description.tsv"
info_path="data/scRNA_metadata.tsv"
#Step locks
step1_lock="keep/op1_seuObj_s1.rda";S1=F;
step2_lock="keep/op1_seuObj_s2.rda";S2=F;
step2_ann="keep/op1_clu_annotation.rda";
#output
output_dir="output/op1"
out_plot_path1="./output/op1/op1_umap_ann.svg"
out_plot_path2="./output/op1/op1_umap_PIEZO1.svg"
out_plot_path3="./output/op1/op1_exp_PIEZO1.svg"
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}else{
  print("Dir already exists!")
}

#Load data ----
geno=fread(geno_path,data.table = F)
see(geno);rownames(geno)=geno[,1];geno=geno[,-1];
header=fread(header_path,data.table = F)
info=fread(info_path,data.table = F);rownames(info)=info$sampleID;
info$subIDmXcond=paste(info$subIDm,info$batchCond,sep="_")

#Check sample order in geno and info ---
ord=seq(1,nrow(info))
matches=match(colnames(geno),info$sampleID)
cat("sample order in geno and info is identical:",all.equal(matches,ord),"\n")

colors=c("red","brown","green","orange","blue","grey")
ann_df=info[,seq(11,18)]
ann_v=rep("ZZZ",nrow(ann_df))

for(k in 1:ncol(ann_df)){
  run_v=ann_df[,k]
  repl_bol=run_v!="ZZZ"
  ann_v[repl_bol]=run_v[repl_bol]
}

### Step 1 ----
if(S1 & S2){
  seuObj <- CreateSeuratObject(counts = geno, 
                             project = "FER", 
                             min.cells = 3, 
                             min.features = 200,
                             meta.data = info)
  
  #Compute percentage of mt
  seuObj[["percent.mt"]] <- PercentageFeatureSet(seuObj, pattern = "^MT-")
  
  #Quality control and filtering of cells based on counts, features and mt
  VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  plot1 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  
  #Filtering
  seuObj <- subset(seuObj, subset = percent.mt < 5)
  
  #Normalization
  seuObj <- NormalizeData(seuObj)
  
  #Feature Selection
  seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)
  
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(seuObj), 10)
  
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(seuObj)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 + plot2
  
  #Scaling
  all.genes <- rownames(seuObj)
  seuObj <- ScaleData(seuObj, features = all.genes)
  
  #PCA
  seuObj <- RunPCA(seuObj, features = VariableFeatures(object = seuObj))
  #Clustering
  seuObj <- JackStraw(seuObj, num.replicate = 100)
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
  n_best_dims=n_pca
  seuObj <- FindNeighbors(seuObj, dims = 1:n_best_dims)
  #I decided to minimize small difficult to annotate clusters
  seuObj <- FindClusters(seuObj, resolution = 0.32)
  head(Idents(seuObj), 5)
  seuObj <- RunUMAP(seuObj, dims = 1:n_best_dims)
  #Replacment of the umap components to the ones of the publicated article
  tmp_umap=seuObj@reductions[["umap"]]@cell.embeddings
  new_umap=cbind(seuObj@meta.data$UMAP1_ALL,seuObj@meta.data$UMAP2_ALL)
  colnames(new_umap)=colnames(tmp_umap)
  rownames(new_umap)=rownames(tmp_umap)
  seuObj@reductions[["umap"]]@cell.embeddings=new_umap
  png("./output/op1_umap.png", units="in", width=22, height=12, res=300)
    DimPlot(seuObj, reduction = "umap", label=T, 
            pt.size = 0.8, label.size = 8)
  dev.off()
  #Preserve this final big step to replicate in the future the same analysis
  save(seuObj,file=step2_lock)
  #Compute annotation of the clusters
  clu_markers <- findmarkergenes(object=seuObj,
                                 species = 'Human',
                                 match_CellMatch = T,
                                 tissue = 'Brain')
  
  clu_ann <- scCATCH(object=clu_markers$clu_markers,
                     species = "Human",
                     cancer = NULL,
                     tissue = "Brain")
  #Save annotation
  save(clu_ann,file=step2_ann)
  
}else{
  load(step2_lock)
  load(step2_ann)
}

#Update the annotation of the clusters based on scCATCH result
orig_clusters=Idents(object = seuObj)
lvs=levels(orig_clusters);names(lvs)=lvs

#Manual overide after having checked for correctess
clu_ann$cell_type[2]=clu_ann$cell_type[1]
clu_ann$cell_type[5]="OPC"
clu_ann$cell_type[6]="Neuron"

lvs[clu_ann[,1]]=clu_ann$cell_type
lvs=paste(names(lvs),lvs,sep=":")
lvs[9]="8:Unidentified"
lvs[10]="9:Unidentified"
names(lvs) <- levels(seuObj)
seuObj <- RenameIdents(seuObj, lvs)

#PLots
svg(filename=out_plot_path1, width=15, height=15, pointsize=12)
DimPlot(seuObj, reduction = "umap", label=T, pt.size = 0.8, label.size = 8)
dev.off()

svg(out_plot_path2, width=15, height=15, pointsize=12)
  FeaturePlot(seuObj, features = c("PIEZO1"), pt.size = 2.5,
              min.cutoff = 1, max.cutoff = 3, cols=c("grey85","darkorange1"),
              order=T)
dev.off()

# svg(out_plot_path3, width=15, height=15, pointsize=12)
# VlnPlot(seuObj, features = c("PIEZO1"), split.by = "batchCond", split.plot = TRUE,
#         cols=c("yellow","grey"))
# dev.off()

#Data fetching ----
batchesXP = FetchData(object = seuObj, vars = c('batchCond', 'ident','PIEZO1'))
save(batchesXP,file="./output/op1_batchesXP.rda")

subidXP = FetchData(object = seuObj, vars = c('subIDmXcond', 'ident','PIEZO1'))
save(subidXP,file="./output/op1_subidXP.rda")

geno = seuObj@assays[["RNA"]]
info = seuObj@meta.data
info$cls = Idents(object = seuObj)
save(geno,info,file="./output/op1_data_ready.rda")















