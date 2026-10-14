rm(list = ls())
library(Seurat)
library(Matrix)
library(tidyverse)
library(org.Hs.eg.db)
library(DoubletFinder)
library(GEOquery)

smps = list.files("rawdata/GSE214695/", "mtx.gz")
smps = str_remove(smps, "_matrix.mtx.gz")

sc.obj.list = list()
for(s in smps){
  print(s)
  sc = ReadMtx(mtx = paste0("rawdata/GSE214695/",s,"_matrix.mtx.gz"),
               cells = paste0("rawdata/GSE214695/",s,"_barcodes.tsv.gz"),
               features = paste0("rawdata/GSE214695/",s,"_features.tsv.gz"))
  s = str_match(s, "GSM[0-9]*")[1]
  sc = CreateSeuratObject(sc, project = s)
  sc = subset(sc, nCount_RNA > 0)
  sc.obj.list[[s]] = sc
}

sc.obj.list = lapply(sc.obj.list, function(x){
  x[["percent.mt"]] = PercentageFeatureSet(x, pattern = "^MT-")
  return(x)
})

sc.obj.list = lapply(sc.obj.list, function(x){
  x = subset(x, nFeature_RNA > 200 & percent.mt < 25)
  return(x)
})

for(i in names(sc.obj.list)){
  WGCNA::collectGarbage()
  x = sc.obj.list[[i]]
  x = NormalizeData(x)
  x = FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x = ScaleData(x)
  x = RunPCA(x)
  x = RunUMAP(x, dims = 1:10)
  x = FindNeighbors(x, reduction = "pca", dims = 1:10)
  x = FindClusters(x, resolution = 0.5)
  
  sweep.res.list = paramSweep(x, PCs = 1:10)
  sweep.stats = summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn = find.pK(sweep.stats)
  
  homotypic.prop = modelHomotypic(x$seurat_clusters)
  nExp_poi = round(0.05*nrow(x@meta.data))
  nExp_poi.adj = round(nExp_poi*(1-homotypic.prop))
  reuse.pann = colnames(x@meta.data) %>% dplyr::last()
  
  x = doubletFinder(x, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, 
                    reuse.pANN = FALSE, sct = FALSE)
  x = doubletFinder(x, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, 
                    reuse.pANN = reuse.pann, sct = FALSE)
  res.doublt = x@meta.data[ncol(x@meta.data)]
  x$"DF.class" = res.doublt
  x = subset(x, DF.class == "Singlet")
  sc.obj.list[[i]] = x
}

sc.obj.merge = merge(sc.obj.list[[1]], sc.obj.list[-1])
sc.obj.merge = NormalizeData(sc.obj.merge)
sc.obj.merge = FindVariableFeatures(sc.obj.merge)
sc.obj.merge = ScaleData(sc.obj.merge)
sc.obj.merge = RunPCA(sc.obj.merge)

sc.obj.merge = IntegrateLayers(object = sc.obj.merge, method = CCAIntegration, 
                               orig.reduction = "pca", 
                               new.reduction = "integrated.cca",
                               verbose = TRUE)

sc.obj.merge[["RNA"]] = JoinLayers(sc.obj.merge[["RNA"]])

df.meta = getGEO(file = "rawdata/GSE214695_series_matrix.txt.gz", getGPL = F)
df.meta = pData(df.meta)
df.meta$group = c("Crohn’s disease" = "CD",
                  "healthy control" = "NC",
                  "ulcerative colitis" = "UC")[df.meta$`disease state:ch1`]
colnames(df.meta)[2] = "sample"
df.meta = df.meta[c("sample","group")]
rownames(df.meta) = df.meta$sample

metadata = sc.obj.merge@meta.data
metadata$Group = df.meta[metadata$orig.ident,]$group
metadata$Sample = metadata$orig.ident
sc.obj.merge = AddMetaData(sc.obj.merge, metadata)
base::saveRDS(sc.obj.merge, "interm/integrated.rds")
