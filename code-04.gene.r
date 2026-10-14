rm(list = ls())
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggsci)
library(ggpubr)
library(msigdbr)
library(IOBR)

gsea.plot = function(res.kegg, top.hall, gene){
  gsdata <- do.call(rbind, lapply(top.hall, enrichplot:::gsInfo, object = res.kegg))
  gsdata$Description = factor(gsdata$Description, levels = top.hall)
  p1 = ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(14) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_manual(values = c(pal_aaas()(10),"gold4")) +
    #scale_color_aaas() +
    ggtitle(gene) +
    geom_hline(yintercept = 0, color = "black", size = 0.8) +
    geom_line(aes_(y = ~runningScore, color = ~Description), size = 1) +
    theme(legend.position = "right", legend.title = element_blank(), legend.background = element_rect(fill = "transparent")) +
    ylab("Running Enrichment Score") + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.line.x = element_blank())
  i = 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1 }
  p2 = ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) + 
    theme_classic(14) + theme(legend.position = "none", 
                              axis.ticks = element_blank(), 
                              axis.text = element_blank(), 
                              axis.line.x = element_blank()) + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    # scale_color_aaas() +
    scale_color_manual(values = c(pal_aaas()(10),"gold4"))
  p = aplot::insert_bottom(p1, p2, height = 0.15)
  return(p)
}

df.m = msigdbr()
df.kegg = subset(df.m, gs_subcollection %in% c("CP:KEGG_LEGACY"))[c(9,1)]

df.exp = read.delim2("cleandata/GSE87466.exp.xls", row.names = 1) |> 
  type.convert(as.is = T)
df.exp = t(df.exp) |> as.data.frame()

gs = read.delim2("02.Target/06.Intersection.xls")
gs = gs$gene

for(g in gs){
  df.sub = df.exp[-which(colnames(df.exp)==g)]
  exp.g = df.exp[[g]]
  res.cor = cor(df.sub, exp.g, method = "spearman") |> as.vector()
  names(res.cor) = colnames(df.sub)
  res.cor = res.cor[order(res.cor, decreasing = T)]
  res.cor = na.omit(res.cor)
  
  res.gsea = GSEA(res.cor, TERM2GENE = df.kegg, pvalueCutoff = 0.05, 
                  seed = 1, pAdjustMethod = "BH", eps = 0)
  top.kegg = res.gsea@result
  write.table(top.kegg, paste0("04.Gene/01.GSEA/", g, ".xls"), quote = F, sep = "\t", row.names = F)
  top.kegg = subset(top.kegg, p.adjust<0.05)
  top.kegg = top.kegg[order(top.kegg$p.adjust, decreasing = F),]
  top.kegg = rbind(subset(top.kegg, NES > 0)[1:5,],
                   subset(top.kegg, NES < 0)[1:5,])
  top.kegg = top.kegg[order(top.kegg$NES, decreasing = T),]
  top.kegg = na.omit(top.kegg)
  top.kegg = top.kegg$Description
  
  p = gsea.plot(res.gsea, top.kegg, g)
  fn1 = paste0("04.Gene/01.GSEA/", g, ".png")
  fn2 = paste0("04.Gene/01.GSEA/", g, ".pdf")
  ggsave(fn1, p, width = 12, height = 6, units = "in", dpi = 300, bg = "white")
  ggsave(fn2, p, width = 12, height = 6, units = "in", dpi = 300, bg = "white")
}



rm(list = ls())
library(IOBR)
df.exp = read.delim2("cleandata/GSE87466.exp.xls", row.names = 1) |> 
  type.convert(as.is = T)
df.cb = deconvo_cibersort(df.exp, arrays = T, absolute = F, perm = 1)
df.cb = df.cb[1:23]
colnames(df.cb) = str_remove(colnames(df.cb), "_CIBERSORT")
colnames(df.cb)[1] = "sample"
write.table(df.cb, "interm/GSE87466.cibersort.xls", quote = F, sep = "\t", row.names = F)



rm(list = ls())
df.im = read.delim2("interm/GSE87466.cibersort.xls", row.names = 1) |> 
  type.convert(as.is = T)
df.grp = read.delim2("cleandata/GSE87466.pheno.xls")
df.im = df.im[df.grp$sample,]

p = apply(df.im, 2, \(x){wilcox.test(x~df.grp$group)$p.value})
df.diff = data.frame(cell = colnames(df.im),
                     p.val = p,
                     p.adj = p.adjust(p, method = "BH"))
df.diff = df.diff[order(df.diff$p.adj, decreasing = F),]
write.table(df.diff, "04.Gene/02.Immune.xls", quote = F, sep = "\t", row.names = F)


df.diff = na.omit(df.diff)
df.diff = df.diff[df.diff$p.adj < 0.05,]
df.im = df.im[df.diff$cell]
df.im$sample = rownames(df.im)
df.im = reshape2::melt(df.im, id.var = "sample", variable.name = "cell")
df.im$group = df.grp$group[match(df.im$sample, df.grp$sample)]
df.diff$lab = paste0(df.diff$cell, "\np.adj = ", signif(df.diff$p.adj,3))
df.im$lab = df.diff$lab[match(df.im$cell, df.diff$cell)]


ggplot(df.im, aes(x = group, y = value)) + 
  geom_violin(aes(fill = group), alpha = 0.5) +
  geom_boxplot(aes(fill = group), width = 0.1) + 
  facet_wrap(~lab, scale = "free", ncol = 3) +
  scale_fill_bmj() + 
  theme_classic2() + 
  labs(x = NULL, y = "Score") +
  theme(axis.text.x.bottom = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        strip.background = element_rect(color = NA),
        legend.position = "top",
        legend.title = element_blank())
ggsave("04.Gene/03.Immune.png", width = 8, height = 11, units = "in", dpi = 300, bg = "white")
ggsave("04.Gene/03.Immune.pdf", width = 8, height = 11, units = "in", dpi = 300, bg = "white")


for(cell in unique(df.im$cell)){
  d = df.im[df.im$cell == cell,]
  ggplot(d, aes(x = group, y = value)) + 
    geom_violin(aes(fill = group), alpha = 0.5) +
    geom_boxplot(aes(fill = group), width = 0.1) + 
    facet_wrap(~lab, scale = "free", ncol = 3) +
    scale_fill_bmj() + 
    theme_classic2() + 
    labs(x = NULL, y = "Score") +
    theme(axis.text.x.bottom = element_blank(),
          axis.ticks.x.bottom = element_blank(),
          strip.background = element_rect(color = NA),
          legend.position = "top",
          legend.title = element_blank())
  fn1 = paste0("04.Gene/03.Immune.Single/", cell, ".png")
  fn2 = paste0("04.Gene/03.Immune.Single/", cell, ".pdf")
  ggsave(fn1, width = 4, height = 5, units = "in", dpi = 300, bg = "white")
  ggsave(fn2, width = 4, height = 5, units = "in", dpi = 300, bg = "white")
}


rm(list = ls())
df.im = read.delim2("interm/GSE87466.cibersort.xls", row.names = 1) |> 
  type.convert(as.is = T)
df.diff = read.delim2("04.Gene/02.Immune.xls") |> type.convert(as.is = T)
df.diff = df.diff[df.diff$p.adj < 0.05,]
df.im = df.im[df.diff$cell]

df.exp = read.delim2("cleandata/GSE87466.exp.xls", row.names = 1) |> 
  type.convert(as.is = T)
df.exp = t(df.exp) |> as.data.frame()
gs = read.delim2("02.Target/06.Intersection.xls")
df.exp = df.exp[gs$gene]
df.exp = df.exp[rownames(df.im),]

res.cor = linkET::correlate(df.exp, df.im, method = "spearman")
df.r = res.cor$r |> as.data.frame()
df.p = res.cor$p |> as.data.frame()
df.r$gene = rownames(df.r)
df.p$gene = rownames(df.p)
df.r = reshape2::melt(df.r, id.var = "gene", variable.name = "cell", value.name = "r")
df.p = reshape2::melt(df.p, id.var = "gene", variable.name = "cell", value.name = "p")
df.cor = merge(df.r, df.p)
write.table(df.cor, "04.Gene/04.Immune.Cor.xls", quote = F, sep = "\t", row.names = F)

df.cor = df.cor[df.cor$p < 0.05,]
df.cor$p[df.cor$p==0] = 1e-14
ggplot(df.cor, aes(x = gene, y = cell)) + 
  geom_point(aes(size = p, color = r)) + 
  scale_size_continuous(range = c(7,3), transform = "log10") + 
  scale_color_gradient2(low = "darkblue", high = "darkred") + 
  labs(x = "Hub genes", y = "Immune cells") +
  SCP::theme_scp()
ggsave("04.Gene/05.Immune.Cor.png", width = 7, height = 8, units = "in", dpi = 300, bg = "white")
ggsave("04.Gene/05.Immune.Cor.pdf", width = 7, height = 8, units = "in", dpi = 300, bg = "white")



rm(list = ls())
df.exp = read.delim2("cleandata/GSE87466.exp.xls", row.names = 1) |> 
  type.convert(as.is = T)
df.exp = t(df.exp) |> as.data.frame()
gs = read.delim2("02.Target/06.Intersection.xls")
df.hub = df.exp[gs$gene]
df.im = df.exp[c("IL1B","IL6","IL10","IL15","IL27","IL33","IFNG","TNF")]

df.grp = read.delim2("cleandata/GSE87466.pheno.xls")
df.im = df.im[df.grp$sample,]
df.diff = read.delim2("02.Target/03.DEG.xls") |> type.convert(as.is = T)
df.diff = df.diff[df.diff$Gene %in% colnames(df.im),]

df.im = df.im[df.diff$Gene]
df.im$sample = rownames(df.im)
df.im = reshape2::melt(df.im, id.var = "sample", variable.name = "gene")
df.im$group = df.grp$group[match(df.im$sample, df.grp$sample)]
df.diff$lab = paste0(df.diff$Gene, "\np.adj = ", signif(df.diff$adj.P.Val,3))
df.im$lab = df.diff$lab[match(df.im$gene, df.diff$Gene)]


ggplot(df.im, aes(x = group, y = value)) + 
  geom_violin(aes(fill = group), alpha = 0.5) +
  geom_boxplot(aes(fill = group), width = 0.1) + 
  facet_wrap(~lab, scale = "free", ncol = 2) +
  scale_fill_bmj() + 
  theme_classic2() + 
  labs(x = NULL, y = "Score") +
  theme(axis.text.x.bottom = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        strip.background = element_rect(color = NA),
        legend.position = "top",
        legend.title = element_blank())
ggsave("04.Gene/06.Cytokines.png", width = 5, height = 8, units = "in", dpi = 300, bg = "white")
ggsave("04.Gene/06.Cytokines.pdf", width = 5, height = 8, units = "in", dpi = 300, bg = "white")


for(gene in unique(df.im$gene)){
  d = df.im[df.im$gene == gene,]
  ggplot(d, aes(x = group, y = value)) + 
    geom_violin(aes(fill = group), alpha = 0.5) +
    geom_boxplot(aes(fill = group), width = 0.1) + 
    facet_wrap(~lab, scale = "free", ncol = 3) +
    scale_fill_bmj() + 
    theme_classic2() + 
    labs(x = NULL, y = "Expression") +
    theme(axis.text.x.bottom = element_blank(),
          axis.ticks.x.bottom = element_blank(),
          strip.background = element_rect(color = NA),
          legend.position = "top",
          legend.title = element_blank())
  fn1 = paste0("04.Gene/06.Cytokines.Single/", gene, ".png")
  fn2 = paste0("04.Gene/06.Cytokines.Single/", gene, ".pdf")
  ggsave(fn1, width = 4, height = 5, units = "in", dpi = 300, bg = "white")
  ggsave(fn2, width = 4, height = 5, units = "in", dpi = 300, bg = "white")
}










df.im = df.exp[df.diff$Gene]
res.cor = linkET::correlate(df.hub, df.im, method = "spearman")
df.r = res.cor$r |> as.data.frame()
df.p = res.cor$p |> as.data.frame()
df.r$hub = rownames(df.r)
df.p$hub = rownames(df.p)
df.r = reshape2::melt(df.r, id.var = "hub", variable.name = "cytokine", value.name = "r")
df.p = reshape2::melt(df.p, id.var = "hub", variable.name = "cytokine", value.name = "p")
df.cor = merge(df.r, df.p)
write.table(df.cor, "04.Gene/07.Cytokine.Cor.xls", quote = F, sep = "\t", row.names = F)


df.cor = df.cor[df.cor$p < 0.05,]
df.cor$p[df.cor$p==0] = 1e-14
ggplot(df.cor, aes(x = hub, y = cytokine)) + 
  geom_point(aes(size = p, color = r)) + 
  scale_size_continuous(range = c(8,4), transform = "log10") + 
  scale_color_gradient2(low = "darkblue", high = "darkred") + 
  labs(x = "Hub genes", y = "Cytokines") +
  SCP::theme_scp()
ggsave("04.Gene/08.Cytokine.Cor.png", width = 6, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("04.Gene/08.Cytokine.Cor.pdf", width = 6, height = 5, units = "in", dpi = 300, bg = "white")


rm(list = ls())
df.mir = read.csv("04.Gene/11.miRNA.csv")
df.mir = df.mir[1:3]
df.mir = df.mir[!duplicated(df.mir),]
df.cnt = table(df.mir$ID) |> data.frame()
df.cnt = df.cnt[df.cnt$Freq == 4,]
cat(as.character(df.cnt$Var1), sep = "\n")
df.mir = df.mir[df.mir$ID %in% df.cnt$Var1,]
df.mir = df.mir[c(1,3)]

df.lnc = read.csv("04.Gene/12.lncRNA.csv")
df.lnc = df.lnc[c(1,3)]
df.lnc = df.lnc[!duplicated(df.lnc),]
df.cnt = table(df.lnc$Target) |> data.frame()
df.cnt = df.cnt[order(df.cnt$Freq, decreasing = T),]
df.cnt = df.cnt[1:10,]
df.lnc = df.lnc[df.lnc$Target %in% df.cnt$Var1,]

colnames(df.lnc) = colnames(df.mir) = c("source","target")
df.full = rbind(df.lnc, df.mir)
write.table(df.full, "interm/cerna.tsv", quote = F, sep = "\t", row.names = F)
write.table(df.lnc['target'], "interm/lncs.tsv", quote = F, sep = "\t", row.names = F)






