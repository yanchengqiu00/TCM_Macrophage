rm(list = ls())
library(Seurat)
library(ggsci)
library(tidyverse)
library(AUCell)
library(gghalves)
library(ggpubr)

sc.obj.merge = readRDS("interm/integrated.rds")
p1 = VlnPlot(sc.obj.merge, "nFeature_RNA", group.by = "orig.ident") + 
  NoLegend() + 
  scale_fill_d3() +
  xlab(NULL) +
  theme(axis.text.x.bottom = element_blank())
p2 = VlnPlot(sc.obj.merge, "nCount_RNA", group.by = "orig.ident") + 
  NoLegend() + 
  scale_fill_d3() +
  xlab(NULL) +
  theme(axis.text.x.bottom = element_blank())
p3 = VlnPlot(sc.obj.merge, "percent.mt", group.by = "orig.ident") + 
  scale_fill_d3() +
  NoLegend()
p = p1+p2+p3
ggsave("01.SC/01.QC.png", p, width = 9, height = 9, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/01.QC.pdf", p, width = 9, height = 9, units = "in", dpi = 300, bg = "white")


dims = 1:50
sc.obj.merge = FindNeighbors(sc.obj.merge, reduction = "integrated.cca", dims = dims)
sc.obj.merge = FindClusters(sc.obj.merge, resolution = 2, algorithm = 4)

sc.obj.merge = RunUMAP(sc.obj.merge, dims = 1:50, reduction = "integrated.cca",
                       n.neighbors = 8, min.dist = 0.5)
DimPlot(sc.obj.merge, label = F, repel = T, shuffle = T)
ggsave("01.SC/02.Cluster.png", width = 6.7, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/02.Cluster.pdf", width = 6.7, height = 6, units = "in", dpi = 300, bg = "white")

DotPlot(sc.obj.merge, cols = c("grey100","darkred"), scale.max = 80, scale.min = 0,
        scale.by = "size", dot.scale = 6, col.min = 0,
        features = c("TOP2A", "PCNA", #TA
                     "SLC26A2", "FABP1", # colonocyte
                     "TRPM5", "SH2D6", # Tuft
                     "MUC2", "TFF3", # goblet
                     "DCN","COL1A2",
                     "PECAM1","CLDN5", # endothelial
                     "S100B", "CDH19" , #glial
                     "CD3D","CD3E",
                     "NKG7","KLRD1",
                     "FOXP3","IL2RA",
                     "CD79A","MS4A1",
                     "CTSW", "KLRB1", # mait
                     "CD38", "MZB1",
                     "CD68","C1QA", # macro
                     "FCGR3B", "CSF3R", # neutrophil
                     "MS4A2","TPSAB1")) +
  coord_flip()
ggsave("01.SC/03.Cluster.Gene.png", width = 11, height = 8, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/03.Cluster.Gene.pdf", width = 11, height = 8, units = "in", dpi = 300, bg = "white")

clst.type = rep("", length = nlevels(sc.obj.merge$seurat_clusters))
names(clst.type) = seq_len(nlevels(sc.obj.merge$seurat_clusters))-1
clst.type[c("0","1","2","4","13","17","20","24")] = "Plasma"
clst.type[c("5")] = "B"
clst.type[c("3","8","12","28","25")] = "T"
clst.type[c("9")] = "Treg"
clst.type[c("6","7","22")] = "NKT"
clst.type[c("33")] = "MAIT"
clst.type[c("18","23")] = "Macrophage"
clst.type[c("29")] = "Neutrophil"
clst.type[c("27")] = "Mast"
clst.type[c("16","19","26")] = "Fibroblast"
clst.type[c("31")] = "Endothelial"
clst.type[c("14")] = "TA"
clst.type[c("21","32")] = "Goblet"
clst.type[c("10","11","15")] = "Colonocyte"
clst.type[c("30")] = "Tuft"
clst.type[c("34")] = "Glial"
which(clst.type == "")

cell.type = sc.obj.merge@meta.data['seurat_clusters']
cell.type$type = clst.type[cell.type$seurat_clusters]
cells = c("Colonocyte","Tuft","Glial","TA","Goblet",
          "Endothelial","Fibroblast","Macrophage",
          "Neutrophil","Mast",
          "Plasma","B","T","Treg","NKT","MAIT")
cell.type$type = factor(cell.type$type, levels = cells)

sc.obj.merge = AddMetaData(sc.obj.merge, metadata = cell.type$type, 
                           col.name = "cell.type")
Idents(sc.obj.merge) = sc.obj.merge$cell.type

cell.cols = rev(pal_d3("category20c")(16))
DimPlot(sc.obj.merge, cols = cell.cols, reduction = "umap", label = F) + 
  theme(aspect.ratio = 1)
ggsave("01.SC/04.Cell.Type.png", width = 6.5, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/04.Cell.Type.pdf", width = 6.5, height = 6, units = "in", dpi = 300, bg = "white")

DotPlot(sc.obj.merge, cols = c("grey100","darkred"), scale.max = 80, scale.min = 0,
        scale.by = "size", dot.scale = 6, col.min = 0,
        features = c("SLC26A2", "FABP1",
                     "TRPM5", "SH2D6",
                     "S100B", "CDH19",
                     "TOP2A", "PCNA", 
                     "MUC2", "TFF3",
                     "PECAM1","CLDN5",
                     "DCN","COL1A2",
                     "CD68","C1QA", 
                     "FCGR3B","CSF3R", 
                     "MS4A2","TPSAB1",
                     "CD38", "MZB1",
                     "CD79A","MS4A1",
                     "CD3D","CD3E",
                     "FOXP3","IL2RA",
                     "NKG7","KLRD1",
                     "CTSW", "KLRB1"
                     )) +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))
ggsave("01.SC/05.Annot.Gene.png", width = 12, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/05.Annot.Gene.pdf", width = 12, height = 6, units = "in", dpi = 300, bg = "white")

saveRDS(sc.obj.merge, "interm/annotated.rds")



rm(list = ls())
library(SCP)
library(BiocParallel)
register(SnowParam(workers = 1, progressbar = TRUE))
sc.obj.merge = readRDS("interm/annotated.rds")
sc.obj.merge = RunDEtest(srt = sc.obj.merge, group_by = "cell.type", 
                         fc.threshold = 1, only.pos = FALSE)

DEGs = sc.obj.merge@tools$DEtest_cell.type$AllMarkers_wilcox
DEGs = subset(DEGs, avg_log2FC > 2 & p_val_adj < 0.05 & pct.1 > 0.5)
df.out = DEGs[c(7,6,3,4,1,5)]
colnames(df.out)[1] = "cell"

sc.obj.merge[["RNA"]] = as(sc.obj.merge[["RNA"]], "Assay")
sc.obj.merge = AnnotateFeatures(sc.obj.merge, 
                                species = "Homo_sapiens", 
                                db = c("TF", "CSPA"))

ht = FeatureHeatmap(
  srt = sc.obj.merge, group.by = "cell.type", features = DEGs$gene, 
  feature_split = DEGs$group1,
  cell_annotation_palette = "simspec",
  species = "Homo_sapiens", 
  db = c("KEGG"), 
  anno_terms = TRUE,
  feature_annotation = c("TF", "CSPA"), 
  feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  height = 10, width = 8
)

png("01.SC/06.Cell.Marker.Heatmap.png", width = 18, height = 12, units = "in", bg = "white", res = 300)
print(ht$plot)
dev.off()
pdf("01.SC/06.Cell.Marker.Heatmap.pdf", width = 18, height = 12)
print(ht$plot)
dev.off()




rm(list = ls())
library(scMetabolism)
sc = readRDS("interm/annotated.rds")

res.metab = list()
for(smp in unique(sc$orig.ident)){
  sc.sub = subset(sc, orig.ident == smp)
  countexp = sc.sub@assays$RNA$counts
  rm(sc.sub)
  WGCNA::collectGarbage()
  res.metab[[smp]] = sc.metabolism(as.matrix(countexp), method = "AUCell", ncores = 1)
}
res.metab = purrr::reduce(res.metab, cbind)
df.out = cbind(pwy = rownames(res.metab), res.metab)
saveRDS(df.out, "interm/scmetabolism.rds")



rm(list = ls())
library(pheatmap)
sc = readRDS("interm/annotated.rds")
df.meta = sc@meta.data
df.score = readRDS("interm/scmetabolism.rds")
df.score = df.score[make.names(rownames(df.meta))]

df.score.grp = t(df.score) |> 
  aggregate(by = list(group = df.meta$Group), FUN = mean) |> as.data.frame()
rownames(df.score.grp) = df.score.grp[[1]]
df.score.grp = df.score.grp[-1]

p = pheatmap(t(df.score.grp), scale = "row", cluster_cols = F)
p = ggplotify::as.ggplot(p)
ggsave("01.SC/07.scMetabolism.Group.png", p, width = 7, height = 11, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/07.scMetabolism.Group.pdf", p, width = 7, height = 11, units = "in", dpi = 300, bg = "white")


df.score.cell = t(df.score) |> 
  aggregate(by = list(cell = as.character(df.meta$cell.type)), FUN = mean) |> as.data.frame()
rownames(df.score.cell) = df.score.cell[[1]]
df.score.cell = df.score.cell[-1]
p = pheatmap(t(df.score.cell), scale = "row")
p = ggplotify::as.ggplot(p)
ggsave("01.SC/08.scMetabolism.Cell.png", p, width = 9, height = 11, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/08.scMetabolism.Cell.pdf", p, width = 9, height = 11, units = "in", dpi = 300, bg = "white")


df.score = t(df.score) |> as.data.frame()
df.score = df.score['Riboflavin metabolism']
df.score$cell = rownames(df.score)
rownames(df.meta) = make.names(rownames(df.meta))
df.score$group = df.meta[df.score$cell,'Group']
df.score$cell.type = df.meta[df.score$cell, 'cell.type']
colnames(df.score)[1] = "value"
cells = c("Colonocyte","Tuft","Glial","TA","Goblet",
          "Endothelial","Fibroblast","Macrophage",
          "Neutrophil","Mast",
          "Plasma","B","T","Treg","NKT","MAIT")
df.score$cell.type = factor(df.score$cell.type, levels = cells)

ggplot(df.score, aes(x = cell.type, y = value)) + 
  geom_boxplot(aes(color = cell.type)) +
  stat_compare_means(label.x = 8.5, hjust = 0.5) +
  scale_color_manual(values = rev(pal_d3("category20c")(16))) + 
  SCP::theme_scp() +
  theme(legend.position = "none", axis.text.x = element_text(hjust = 1, angle = 30)) + 
  labs(x = "Cell", y = "Riboflavin metabolism")
ggsave("01.SC/09.Selected.Cell.png", width = 9, height = 4.5, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/09.Selected.Cell.pdf", width = 9, height = 4.5, units = "in", dpi = 300, bg = "white")

df.score = df.score[df.score$cell.type == "Macrophage",]
ggplot(df.score, aes(x = group, y = value)) + 
  geom_violin(aes(fill = group), alpha = 0.5) +
  geom_boxplot(aes(fill = group), width = 0.1) + 
  stat_compare_means(label.x = 1.5, hjust = 0.5) +
  scale_fill_bmj() + 
  SCP::theme_scp() +
  theme(legend.position = "none", axis.text.x = element_text(hjust = 1, angle = 30)) + 
  labs(x = "Group", y = "Riboflavin metabolism")
ggsave("01.SC/10.Selected.Group.png", width = 3.5, height = 4.5, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/10.Selected.Group.pdf", width = 3.5, height = 4.5, units = "in", dpi = 300, bg = "white")



rm(list = ls())
library(hdWGCNA)
sc.obj.merge = readRDS("interm/annotated.rds")
sc.obj.merge = subset(sc.obj.merge, cell.type == "Macrophage")

sc.obj.merge = SetupForWGCNA(sc.obj.merge, 
                             gene_select = "fraction",
                             fraction = 0.05,
                             wgcna_name = "wgcna")

sc.obj.merge = MetacellsByGroups(sc.obj.merge, 
                                 group.by = c("orig.ident","cell.type"),
                                 reduction = "integrated.cca",
                                 ident.group = "cell.type")
sc.obj.merge = NormalizeMetacells(sc.obj.merge)
sc.obj.merge = SetDatExpr(sc.obj.merge, 
                          group_name = unique(sc.obj.merge$cell.type),
                          group.by = "cell.type",
                          assay = "RNA",
                          slot = "data")
sc.obj.merge = TestSoftPowers(sc.obj.merge, 
                              powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),
                              networkType = 'signed')

plot_list = PlotSoftPowers(sc.obj.merge)
patchwork::wrap_plots(plot_list, ncol=2)
ggsave("01.SC/11.SoftPower.png", width = 7, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/11.SoftPower.pdf", width = 7, height = 7, units = "in", dpi = 300, bg = "white")

sc.obj.merge = ConstructNetwork(sc.obj.merge, 
                                tom_outdir = "interm/TOM",
                                tom_name = "wgcna",
                                overwrite_tom = T)
sc.obj.merge = ModuleEigengenes(sc.obj.merge, 
                                group.by.vars = "orig.ident")
sc.obj.merge = ModuleConnectivity(sc.obj.merge,
                                  group.by = "cell.type",
                                  group_name = unique(sc.obj.merge$cell.type))

png("01.SC/12.Dendrogram.png", width = 6, height = 4, units = "in", bg = "white", res = 300)
PlotDendrogram(sc.obj.merge, main='hdWGCNA Dendrogram')
dev.off()
pdf("01.SC/12.Dendrogram.pdf", width = 6, height = 4)
PlotDendrogram(sc.obj.merge, main='hdWGCNA Dendrogram')
dev.off()

png("01.SC/13.KME.png", width = 6, height = 7, units = "in", bg = "white", res = 300)
PlotKMEs(sc.obj.merge, ncol = 3)
dev.off()
pdf("01.SC/13.KME.pdf", width = 6, height = 7)
PlotKMEs(sc.obj.merge, ncol = 3)
dev.off()

plot_list = ModuleFeaturePlot(sc.obj.merge, features='hMEs', order=TRUE)
patchwork::wrap_plots(plot_list, ncol=3)
ggsave("01.SC/14.Module.Feature.png", width = 7, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/14.Module.Feature.pdf", width = 7, height = 6, units = "in", dpi = 300, bg = "white")

hMEs = GetMEs(sc.obj.merge)
hMEs = cbind(cell = rownames(hMEs), hMEs)
write.table(hMEs, "interm/hdwgcna.xls", quote = F, sep = "\t", row.names = F)

modules = GetModules(sc.obj.merge) |> subset(module != 'grey')
write.table(modules, "01.SC/15.Gene.Module.xls", quote = F, sep = "\t", row.names = F)



rm(list = ls())
sc = readRDS("interm/annotated.rds")
df.meta = sc@meta.data
df.score = readRDS("interm/scmetabolism.rds")
df.score = df.score[-1]
df.score = t(df.score) |> as.data.frame()

df.wgcna = read.delim2("interm/hdwgcna.xls", row.names = 1) |> type.convert(as.is = T)
rownames(df.wgcna) = make.names(rownames(df.wgcna))
df.score = df.score[rownames(df.wgcna),]
df.score = df.score['Riboflavin metabolism']
res.cor = linkET::correlate(df.wgcna, df.score, method = "spearman")
df.r = res.cor$r |> as.data.frame()
df.p = res.cor$p |> as.data.frame()
df.r = cbind(module = rownames(df.r), df.r)
df.p = cbind(module = rownames(df.p), df.p)
colnames(df.r)[2] = "r"
colnames(df.p)[2] = "p"
df.cor = merge(df.r, df.p)
df.cor = df.cor[df.cor$module != "grey",]
write.table(df.cor, "01.SC/15.Score.Module.xls", quote = F, sep = "\t", row.names = F)

ggplot(df.cor, aes(x = 1, y = module)) + 
  geom_point(aes(color = r), size = 10) + 
  geom_text(aes(label = signif(r,2))) +
  scale_color_gradient2(low = "blue", high = "red") +
  theme_classic2() +
  labs(x = "Riboflavin\nmetabolism", y = "Module") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave("01.SC/17.Score.Module.png", width = 3, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("01.SC/17.Score.Module.pdf", width = 3, height = 6, units = "in", dpi = 300, bg = "white")


