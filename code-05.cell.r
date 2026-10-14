rm(list = ls())
library(tidyverse)
library(SCP)
library(Seurat)
library(harmony)
library(ggsci)
library(monocle3)
library(ClusterGVis)
library(org.Hs.eg.db)
library(CellChat)
library(patchwork)

sc.sub = readRDS("interm/annotated.rds")
sc.sub = subset(sc.sub, cell.type == "Macrophage")
sc.sub = NormalizeData(sc.sub)
sc.sub = FindVariableFeatures(sc.sub)
sc.sub = ScaleData(sc.sub)
sc.sub = RunPCA(sc.sub)
sc.sub = RunHarmony(sc.sub, "orig.ident")
sc.sub = FindNeighbors(sc.sub, dims = 1:30, reduction = "harmony")
sc.sub = FindClusters(sc.sub, resolution = 0.1)

res.marker = FindAllMarkers(sc.sub, only.pos = T)
res.marker = res.marker[c(6,7,2,3,4,1,5)]
res.marker = res.marker[res.marker$p_val_adj < 0.05, ]

df.meta = sc.sub@meta.data
df.meta$subcluster = c("0" = "Macro_TIMP1",
                      "1" = "Macro_TMSB4X",
                      "2" = "Macro_A2M",
                      "3" = "Macro_CCL19")[as.character(df.meta$seurat_clusters)]
sc.sub = AddMetaData(sc.sub, df.meta['subcluster'])
sc.sub$subcluster = factor(sc.sub$subcluster)

sc.sub = RunUMAP(sc.sub, dims = 1:30, reduction = "harmony", 
                 n.neighbors = 50, min.dist = 0.5)
DimPlot(sc.sub, label = T, group.by = "subcluster", repel = T) + NoLegend()
ggsave("05.Macrophage/01.Subcluster.png", width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("05.Macrophage/01.Subcluster.pdf", width = 5, height = 5, units = "in", dpi = 300, bg = "white")


res.marker$cluster = c("0" = "Macro_TIMP1",
                       "1" = "Macro_TMSB4X",
                       "2" = "Macro_A2M",
                       "3" = "Macro_CCL19")[as.character(res.marker$cluster)]
res.marker$cluster = factor(res.marker$cluster)
df.plot = subset(res.marker, pct.1 > 0.5) |> 
  group_by(cluster) |> 
  top_n(n = 3, wt = -p_val_adj)
df.plot = df.plot[order(df.plot$cluster, decreasing = F),]

DotPlot(sc.sub, c(df.plot$gene), group.by = "subcluster") +
  coord_flip() + 
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))
ggsave("05.Macrophage/02.Top.Marker.png", width = 7, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("05.Macrophage/02.Top.Marker.pdf", width = 7, height = 7, units = "in", dpi = 300, bg = "white")


df.sub = sc.sub@meta.data
df.sub$cell = rownames(df.sub)
df.sub = df.sub[c('cell','subcluster')]
write.table(df.sub, "interm/subcluster.tsv", quote = F, sep = "\t", row.names = F)

Idents(sc.sub) = "subcluster"
sc.sub[["RNA"]] = as(sc.sub[["RNA"]], "Assay")
res.marker = res.marker[res.marker$pct.1 > 0.5, ]
sc.sub = AnnotateFeatures(sc.sub, assays = "RNA",
                          species = "Homo_sapiens", 
                          db = c("TF", "CSPA"))

ht = FeatureHeatmap(
  srt = sc.sub, group.by = "subcluster", features = res.marker$gene, 
  feature_split = res.marker$cluster,
  cell_annotation_palette = "simspec",
  species = "Homo_sapiens", 
  db = c("KEGG"), 
  anno_terms = TRUE,
  feature_annotation = c("TF", "CSPA"), 
  feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  height = 8, width = 7
)

png("05.Macrophage/03.Subcluster.Enrich.png", width = 15, height = 10, units = "in", bg = "white", res = 300)
print(ht$plot)
dev.off()
pdf("05.Macrophage/03.Subcluster.Enrich.pdf", width = 15, height = 10)
print(ht$plot)
dev.off()


gs = read.delim2("02.Target/06.Intersection.xls")
p = VlnPlot(sc.sub, gs$gene, split.by = "Group", ncol = 4) +
  scale_fill_bmj() +
  theme(legend.position = "right")
p[[1]] = p[[1]] + scale_fill_bmj()
p[[2]] = p[[2]] + scale_fill_bmj()
p[[3]] = p[[3]] + scale_fill_bmj()
p
ggsave("05.Macrophage/04.Subcluster.Hub.png", width = 9, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("05.Macrophage/04.Subcluster.Hub.pdf", width = 9, height = 4, units = "in", dpi = 300, bg = "white")

saveRDS(sc.sub, "interm/subcluster.rds")



rm(list = ls())
sc.sub = readRDS("interm/subcluster.rds")
gs = read.delim2("02.Target/06.Intersection.xls")
data <- GetAssayData(sc.sub, assay = 'RNA', layer = 'counts')
cell_metadata <- sc.sub@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data), 
                              row.names = rownames(data))
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 10)
cds <- align_cds(cds, alignment_group = "orig.ident")
cds <- reduce_dimension(cds, preprocess_method = "PCA")
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = F)
cds <- order_cells(cds)

p1 = plot_cells(cds, cell_size = 0.5,
                trajectory_graph_segment_size = 0.5,
                color_cells_by = "subcluster",
                label_cell_groups=F,
                label_leaves=FALSE, 
                label_roots = F,
                label_branch_points=FALSE,
                graph_label_size=4) + 
  theme(aspect.ratio = 1) + 
  scale_color_d3()
p2 = plot_cells(cds, cell_size = 0.5,
                trajectory_graph_segment_size = 0.5,
                color_cells_by = "pseudotime",
                label_cell_groups=FALSE,
                label_leaves=FALSE, label_roots = F,
                label_branch_points=FALSE,
                graph_label_size=4) + 
  theme(aspect.ratio = 1)
p1+p2
ggsave("05.Macrophage/05.Pseudotime.png", width = 8, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("05.Macrophage/05.Pseudotime.pdf", width = 8, height = 4, units = "in", dpi = 300, bg = "white")


cell.nc = Cells(subset(sc.sub, Group == "NC"))
cell.uc = Cells(subset(sc.sub, Group == "UC"))

p1 = plot_cells(cds[,cell.nc], cell_size = 0.5,
           trajectory_graph_segment_size = 0.5,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE, label_roots = F,
           label_branch_points=FALSE,
           graph_label_size=4) + 
  ggtitle("NC") + 
  theme(plot.title = element_text(hjust = 0.5))
p2 = plot_cells(cds[,cell.uc], cell_size = 0.5,
           trajectory_graph_segment_size = 0.5,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE, label_roots = F,
           label_branch_points=FALSE,
           graph_label_size=4) + 
  ggtitle("UC") + 
  theme(plot.title = element_text(hjust = 0.5))
p1+p2
ggsave("05.Macrophage/06.Pseudotime.Group.png", width = 8, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("05.Macrophage/06.Pseudotime.Group.pdf", width = 8, height = 4, units = "in", dpi = 300, bg = "white")

df.time = pseudotime(cds) |> data.frame()
sc.sub = AddMetaData(sc.sub, df.time)
df.meta = sc.sub@meta.data
ggplot(df.meta, aes(x = pseudotime.cds.)) + 
  geom_density(aes(fill = Group, color = Group), alpha = 0.5) +
  scale_fill_bmj() + 
  scale_color_bmj() +
  SCP::theme_scp() + 
  labs(x = "Pseudotime", y = "Density")
ggsave("05.Macrophage/07.Pseudotime.Density.png", width = 8, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("05.Macrophage/07.Pseudotime.Density.pdf", width = 8, height = 4, units = "in", dpi = 300, bg = "white")


plot_genes_in_pseudotime(cds[gs$gene,], ncol = 2,
                         color_cells_by="pseudotime",
                         label_by_short_name = F, min_expr= 0.001)
ggsave("05.Macrophage/08.Pseudotime.Hub.png", width = 8, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("05.Macrophage/08.Pseudotime.Hub.pdf", width = 8, height = 4, units = "in", dpi = 300, bg = "white")



genes = graph_test(cds)
top_genes = row.names(subset(genes, q_value == 0 & morans_I > 0.25))
mat = pre_pseudotime_matrix(cds_obj = cds, gene_list = top_genes)
class(mat) = "matrix"

filter.std = ClusterGVis:::filter.std
fix(filter.std)
ck = clusterData(obj = mat, cluster.method = "kmeans", cluster.num = 5)

enrich = enrichCluster(object = ck,
                       OrgDb = org.Hs.eg.db,
                       type = "BP",
                       organism = "hsa",
                       pvalueCutoff = 0.5,
                       topn = 5,
                       seed = 42)

set.seed(42)
markGenes = sample(rownames(mat),50,replace = F)

png("05.Macrophage/09.Pseudotime.Cluster.png", width = 15, height = 10, units = "in", bg = "white", res = 300)
visCluster(object = ck,
           plot.type = "both",
           add.sampleanno = F,
           annoTerm.data = enrich,
           markGenes = markGenes)
dev.off()
pdf("05.Macrophage/09.Pseudotime.Cluster.pdf", width = 15, height = 10)
visCluster(object = ck,
           plot.type = "both",
           add.sampleanno = F,
           annoTerm.data = enrich,
           markGenes = markGenes)
dev.off()



rm(list = ls())
library(ReactomeGSA)
library(gplots)
sc.sub = readRDS("interm/subcluster.rds")
sc.sub$cellgroup = paste0(sc.sub$Group, "_", sc.sub$subcluster)
Idents(sc.sub) = "cellgroup"
gsva_result = analyse_sc_clusters(sc.sub, verbose = TRUE)
pathway_expression = pathways(gsva_result)

png("05.Macrophage/10.ReactomeGSA.png", width = 10, height = 10, units = "in", bg = "white", res = 300)
par(oma = c(3,1,5,1))
plot_gsva_heatmap(gsva_result, max_pathways = 30, margins = c(10,30), 
                  truncate_names = F, col = bluered, lhei = c(2, 10))
dev.off()
pdf("05.Macrophage/10.ReactomeGSA.pdf", width = 12, height = 12)
par(oma = c(3,1,5,1))
plot_gsva_heatmap(gsva_result, max_pathways = 30, margins = c(10,30), 
                  truncate_names = F, col = bluered, lhei = c(2, 10))
dev.off()





rm(list = ls())
sc.obj.merge = readRDS("interm/annotated.rds")
df.sub = read.delim2("interm/subcluster.tsv")
df.meta = sc.obj.merge@meta.data
df.meta$cell = rownames(df.meta)
df.meta = merge(df.meta, df.sub, by = "cell", all = T)
df.meta$subcluster = ifelse(is.na(df.meta$subcluster), as.character(df.meta$cell.type), df.meta$subcluster)
rownames(df.meta) = df.meta$cell
sc.obj.merge = AddMetaData(sc.obj.merge, df.meta['subcluster'])

cellchat = createCellChat(object = sc.obj.merge, 
                          group.by = "subcluster", 
                          assay = "RNA")
CellChatDB = CellChatDB.human
CellChatDB.use = subsetDB(CellChatDB)
cellchat@DB = CellChatDB.use
cellchat = subsetData(cellchat)
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)
cellchat = smoothData(cellchat, adj = PPI.human)
cellchat = computeCommunProb(cellchat, type = "triMean")
cellchat = computeCommunProbPathway(cellchat)
cellchat = aggregateNet(cellchat)

p1 = netVisual_heatmap(cellchat, measure = "count")
p2 = netVisual_heatmap(cellchat, measure = "weight")

ht_opt$legend_gap = unit(10, "mm")
png("05.Macrophage/11.CellChat.HM.png", width = 9, height = 5, units = "in", bg = "white", res = 300)
draw(p1 + p2, newpage = F, 
     column_title_gp = gpar(fontsize = 14),
     auto_adjust = T, heatmap_legend_side = "right")
dev.off()
pdf("05.Macrophage/11.CellChat.HM.pdf", width = 9, height = 5)
draw(p1 + p2, newpage = F, 
     column_title_gp = gpar(fontsize = 14),
     auto_adjust = T, heatmap_legend_side = "right")
dev.off()


sc.obj.merge = subset(sc.obj.merge, subcluster != "Macro_CCL19")
cellchat.c = createCellChat(object = subset(sc.obj.merge, Group == "NC"), 
                            group.by = "subcluster", assay = "RNA")
cellchat.d = createCellChat(object = subset(sc.obj.merge, Group == "UC"), 
                            group.by = "subcluster", assay = "RNA")

cellchat.c@DB = cellchat.d@DB = CellChatDB.use

cellchat.c = subsetData(cellchat.c)
cellchat.c = identifyOverExpressedGenes(cellchat.c)
cellchat.c = identifyOverExpressedInteractions(cellchat.c)
#cellchat.c = smoothData(cellchat.c, adj = PPI.human)
cellchat.c = computeCommunProb(cellchat.c, type = "triMean")
cellchat.c = filterCommunication(cellchat.c, min.cells = 10)
cellchat.c = computeCommunProbPathway(cellchat.c)
cellchat.c = aggregateNet(cellchat.c)

cellchat.d = subsetData(cellchat.d)
cellchat.d = identifyOverExpressedGenes(cellchat.d)
cellchat.d = identifyOverExpressedInteractions(cellchat.d)
#cellchat.d = smoothData(cellchat.d, adj = PPI.human)
cellchat.d = computeCommunProb(cellchat.d, type = "triMean")
cellchat.d = filterCommunication(cellchat.d, min.cells = 10)
cellchat.d = computeCommunProbPathway(cellchat.d)
cellchat.d = aggregateNet(cellchat.d)

object.list <- list(Control = cellchat.c, Disease = cellchat.d)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

for (i in 1:length(object.list)) {
  object.list[[i]] = netAnalysis_computeCentrality(object.list[[i]])
}

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

png("05.Macrophage/12.Interaction.Diff.png", width = 9, height = 5, units = "in", bg = "white", res = 300)
draw(gg1 + gg2, newpage = F, column_title = "Disease vs Control", 
     column_title_gp = gpar(fontsize = 14),
     auto_adjust = T, heatmap_legend_side = "right")
dev.off()
pdf("05.Macrophage/12.Interaction.Diff.pdf", width = 9, height = 5)
draw(gg1 + gg2, newpage = F, column_title = "Disease vs Control", 
     column_title_gp = gpar(fontsize = 14),
     auto_adjust = T, heatmap_legend_side = "right")
dev.off()



for(i in 7:9){
  cell = colnames(cellchat@net$Control$prob)[i]
  p1 = rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = i, 
               targets.use = NULL, stacked = T, do.stat = TRUE, 
               color.use = c("Control" = "#2A6EBBFF","Disease" = "#F0AB00FF")) + 
    ggtitle(paste0(cell, " as target")) + 
    theme(plot.title = element_text(hjust = 0.5))
  p2 = rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, 
               targets.use = i, stacked = T, do.stat = TRUE, 
               color.use = c("Control" = "#2A6EBBFF","Disease" = "#F0AB00FF")) + 
    ggtitle(paste0(cell, " as target")) + 
    theme(plot.title = element_text(hjust = 0.5))
  p = wrap_plots(plots = list(p1,p2)) + plot_layout(guides = "collect")
  fn1 = paste0("05.Macrophage/13.Pathway.Diff/", cell, ".png")
  fn2 = paste0("05.Macrophage/13.Pathway.Diff/", cell, ".pdf")
}


cell = colnames(cellchat@net$Control$prob)[7]
others = c(1:6,8:18)
p = netVisual_bubble(cellchat, sources.use = 7, targets.use = others, 
                     comparison = c(1,2), angle.x = 90,
                     color.text = c("Control" = "#2A6EBBFF","Disease" = "#F0AB00FF"))
fn1 = paste0("05.Macrophage/14.LR.Diff/", cell, "source.png")
fn2 = paste0("05.Macrophage/14.LR.Diff/", cell, "source.pdf")
ggsave(fn1, p, width = 9, height = 15, units = "in", dpi = 300, bg = "white")
ggsave(fn2, p, width = 9, height = 15, units = "in", dpi = 300, bg = "white")


p = netVisual_bubble(cellchat, sources.use = others, targets.use = 7, 
                     comparison = c(1,2), angle.x = 90, 
                     color.text = c("Control" = "#2A6EBBFF","Disease" = "#F0AB00FF"))
fn1 = paste0("05.Macrophage/14.LR.Diff/", cell, "target.png")
fn2 = paste0("05.Macrophage/14.LR.Diff/", cell, "target.pdf")
ggsave(fn1, p, width = 9, height = 20, units = "in", dpi = 300, bg = "white")
ggsave(fn2, p, width = 9, height = 20, units = "in", dpi = 300, bg = "white")



cell = colnames(cellchat@net$Control$prob)[8]
others = c(1:7,9:18)
p = netVisual_bubble(cellchat, sources.use = 8, targets.use = others, 
                     comparison = c(1,2), angle.x = 90,
                     color.text = c("Control" = "#2A6EBBFF","Disease" = "#F0AB00FF"))
fn1 = paste0("05.Macrophage/14.LR.Diff/", cell, "source.png")
fn2 = paste0("05.Macrophage/14.LR.Diff/", cell, "source.pdf")
ggsave(fn1, p, width = 9, height = 15, units = "in", dpi = 300, bg = "white")
ggsave(fn2, p, width = 9, height = 15, units = "in", dpi = 300, bg = "white")


p = netVisual_bubble(cellchat, sources.use = others, targets.use = 8, 
                     comparison = c(1,2), angle.x = 90, 
                     color.text = c("Control" = "#2A6EBBFF","Disease" = "#F0AB00FF"))
fn1 = paste0("05.Macrophage/14.LR.Diff/", cell, "target.png")
fn2 = paste0("05.Macrophage/14.LR.Diff/", cell, "target.pdf")
ggsave(fn1, p, width = 9, height = 20, units = "in", dpi = 300, bg = "white")
ggsave(fn2, p, width = 9, height = 20, units = "in", dpi = 300, bg = "white")



cell = colnames(cellchat@net$Control$prob)[9]
others = c(1:8,10:18)
p = netVisual_bubble(cellchat, sources.use = 9, targets.use = others, 
                     comparison = c(1,2), angle.x = 90,
                     color.text = c("Control" = "#2A6EBBFF","Disease" = "#F0AB00FF"))
fn1 = paste0("05.Macrophage/14.LR.Diff/", cell, "source.png")
fn2 = paste0("05.Macrophage/14.LR.Diff/", cell, "source.pdf")
ggsave(fn1, p, width = 9, height = 7, units = "in", dpi = 300, bg = "white")
ggsave(fn2, p, width = 9, height = 7, units = "in", dpi = 300, bg = "white")


p = netVisual_bubble(cellchat, sources.use = others, targets.use = 9, 
                     comparison = c(1,2), angle.x = 90, 
                     color.text = c("Control" = "#2A6EBBFF","Disease" = "#F0AB00FF"))
fn1 = paste0("05.Macrophage/14.LR.Diff/", cell, "target.png")
fn2 = paste0("05.Macrophage/14.LR.Diff/", cell, "target.pdf")
ggsave(fn1, p, width = 9, height = 7, units = "in", dpi = 300, bg = "white")
ggsave(fn2, p, width = 9, height = 7, units = "in", dpi = 300, bg = "white")




