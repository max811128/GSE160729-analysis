# load libraries
library(Seurat)
library(tidyverse)

LFD_E1<- Read10X(data.dir = "../Downloads/GSE160729_RAW/LFD_E1/")
LFD_E2<- Read10X(data.dir = "../Downloads/GSE160729_RAW/LFD_E2/")
LFD_E3<- Read10X(data.dir = "../Downloads/GSE160729_RAW/LFD_E3/")

HFD_E1<- Read10X(data.dir = "../Downloads/GSE160729_RAW/HFD_E1/")
HFD_E2<- Read10X(data.dir = "../Downloads/GSE160729_RAW/HFD_E2/")
HFD_E3<- Read10X(data.dir = "../Downloads/GSE160729_RAW/HFD_E3/")

# Create Seurat Object and percent.mt
LFD_E1<- CreateSeuratObject(counts = LFD_E1, project = "LFD", min.cells = 3, min.features = 200)
LFD_E1<- PercentageFeatureSet(LFD_E1, pattern = "^mt-", col.name = "percent.mt")
LFD_E2<- CreateSeuratObject(counts = LFD_E2, project = "LFD", min.cells = 3, min.features = 200)
LFD_E2<- PercentageFeatureSet(LFD_E2, pattern = "^mt-", col.name = "percent.mt")
LFD_E3<- CreateSeuratObject(counts = LFD_E3, project = "LFD", min.cells = 3, min.features = 200)
LFD_E3<- PercentageFeatureSet(LFD_E3, pattern = "^mt-", col.name = "percent.mt")

HFD_E1<- CreateSeuratObject(counts = HFD_E1, project = "HFD", min.cells = 3, min.features = 200)
HFD_E1<- PercentageFeatureSet(HFD_E1, pattern = "^mt-", col.name = "percent.mt")
HFD_E2<- CreateSeuratObject(counts = HFD_E2, project = "HFD", min.cells = 3, min.features = 200)
HFD_E2<- PercentageFeatureSet(HFD_E2, pattern = "^mt-", col.name = "percent.mt")
HFD_E3<- CreateSeuratObject(counts = HFD_E3, project = "HFD", min.cells = 3, min.features = 200)
HFD_E3<- PercentageFeatureSet(HFD_E3, pattern = "^mt-", col.name = "percent.mt")

LFD_E1<- subset(LFD_E1, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)
LFD_E2<- subset(LFD_E2, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)
LFD_E3<- subset(LFD_E3, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)
HFD_E1<- subset(HFD_E1, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)
HFD_E2<- subset(HFD_E2, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)
HFD_E3<- subset(HFD_E3, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)

All.list <- list(LFD_E1=LFD_E1, LFD_E2=LFD_E2, LFD_E3=LFD_E3, HFD_E1=HFD_E1, HFD_E2=HFD_E2, HFD_E3=HFD_E3) 

# normalize and identify variable features for each dataset independently
All.list <- lapply(X = All.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})

features <- SelectIntegrationFeatures(object.list = All.list)
All.anchors <- FindIntegrationAnchors(object.list = All.list, anchor.features = features)
All.combined <- IntegrateData(anchorset = All.anchors)

All.combined <- ScaleData(All.combined, verbose = FALSE)
All.combined <- RunPCA(All.combined, npcs = 50, verbose = FALSE)
All.combined <- RunUMAP(All.combined, reduction = "pca", dims = 1:30)
All.combined <- FindNeighbors(All.combined, reduction = "pca", dims = 1:30)
All.combined <- FindClusters(All.combined, resolution = 0.1)
DimPlot(All.combined, reduction = "umap", label = TRUE)

##以組別分開顏色畫Dimplot在同張圖
DimPlot(All.combined, reduction = "umap", group.by = "orig.ident")
##以組別分開畫Dimplot在不同張圖
DimPlot(All.combined, reduction = "umap", split.by = "orig.ident")






