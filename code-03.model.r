rm(list = ls())
library(tidyverse)
library(mlr3verse)
library(ggsci)
library(kernelshap)
library(shapviz)
library(pROC)

df.train = read.delim2("cleandata/GSE87466.exp.xls", row.names = 1) |> 
  type.convert(as.is = T)
df.pheno = read.delim2("cleandata/GSE87466.pheno.xls")
gs = read.delim2("02.Target/06.Intersection.xls")
df.train = df.train[gs$gene,df.pheno$sample] |> t() |> as.data.frame()
df.train$group = df.pheno$group

tsk = as_task_classif(df.train, target = "group", positive = "Disease")
rsmp_cv5 = rsmp("repeated_cv", folds = 5, repeats = 10)
learners = lrns(c("classif.lda",
                  "classif.log_reg",
                  "classif.naive_bayes","classif.qda", 
                  "classif.rpart", "classif.svm", 
                  "classif.xgboost"), predict_type = "prob")
set.seed(25)
design = benchmark_grid(tsk, learners, rsmp_cv5)
bmr = benchmark(design)

lrn_names = c("LDA","Logistic","NB","QDA","RPART","SVM","XGBoost")
png("03.Model/01.BenchMark.AUC.png", width = 8, height = 5, units = "in", bg = "white", res = 300)
autoplot(bmr, measure = msr("classif.auc")) + 
  scale_x_discrete(label = lrn_names)
dev.off()
pdf("03.Model/01.BenchMark.AUC.pdf", width = 8, height = 5)
autoplot(bmr, measure = msr("classif.auc")) + 
  scale_x_discrete(label = lrn_names)
dev.off()

p = autoplot(bmr, type = "roc")
p = p + scale_color_viridis_d("Learner", end = 0.8, aesthetics = c("color", "fill"),
                              labels = lrn_names) + 
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.1), 
        legend.justification = c(1,0))
ggsave("03.Model/02.BenchMark.ROC.png", p, width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("03.Model/02.BenchMark.ROC.pdf", p, width = 5, height = 5, units = "in", dpi = 300, bg = "white")


learner = lrn("classif.naive_bayes", predict_type = "prob")
learner$train(tsk)
shap = kernelshap(learner$model, type = "raw", X = df.train[1:4], seed = 2025)
names(shap$baseline) = c("Disease","Control")
shap_viz = shapviz(shap, which_class = "Disease")

sv_importance(shap_viz, show_numbers = T, number_size = 4) + 
  SCP::theme_scp()
ggsave("03.Model/03.SHAP.bar.png", width = 6, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("03.Model/03.SHAP.bar.pdf", width = 6, height = 4, units = "in", dpi = 300, bg = "white")

sv_importance(shap_viz, kind = "bee") + 
  SCP::theme_scp()
ggsave("03.Model/04.SHAP.beeswarm.png", width = 6, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("03.Model/04.SHAP.beeswarm.pdf", width = 6, height = 4, units = "in", dpi = 300, bg = "white")

sv_force(shap_viz, row_id = 3) +
  ggtitle(paste0("Force plot for ", rownames(df.train)[3]))
ggsave("03.Model/05.SHAP.force.png", width = 6, height = 2.5, units = "in", dpi = 300, bg = "white")
ggsave("03.Model/05.SHAP.force.pdf", width = 6, height = 2.5, units = "in", dpi = 300, bg = "white")


df.test = read.delim2("cleandata/GSE59071.exp.xls", row.names = 1) |> 
  type.convert(as.is = T)
df.pheno = read.delim2("cleandata/GSE59071.pheno.xls")
df.test = df.test[gs$gene,df.pheno$sample]
df.test = t(df.test) |> as.data.frame()
df.pheno$pred = predict(learner$model, df.test, type = "raw")[,'Disease']

roc = roc(predictor = df.pheno$pred, response = df.pheno$group)
auc.lab = paste0("AUC = ", signif(auc(roc),3))
ggroc(roc, legacy.axes = T, linewidth = 1) + 
  annotate("segment", x = 0, xend = 1, y = 0, yend = 1, linetype = 2, linewidth = 0.75) +
  annotate("text", x = 0.6, y = 0.4, label = auc.lab) +
  SCP::theme_scp() + 
  theme(aspect.ratio = 1)
ggsave("03.Model/06.Test.ROC.png", width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("03.Model/06.Test.ROC.pdf", width = 5, height = 5, units = "in", dpi = 300, bg = "white")



rm(list = ls())
df.exp = read.delim2("cleandata/GSE87466.exp.xls", row.names = 1) |> 
  type.convert(as.is = T)
df.pheno = read.delim2("cleandata/GSE87466.pheno.xls")
gs = read.delim2("02.Target/06.Intersection.xls")
df.exp = df.exp[gs$gene,df.pheno$geo_accession]

df.diff = read.delim2("02.Target/03.DEG.xls") |> type.convert(as.is = T)
df.diff$lab = paste0(df.diff$Gene, "\np.adj = ", signif(df.diff$adj.P.Val,3))

df.exp$gene = rownames(df.exp)
df.exp = reshape2::melt(df.exp, id.var = "gene", variable.name = "sample")
df.exp$group = df.pheno$group[match(df.exp$sample, df.pheno$geo_accession)]
df.exp$lab = df.diff$lab[match(df.exp$gene, df.diff$Gene)]

ggplot(df.exp, aes(x = group, y = value)) + 
  geom_boxplot(aes(color = group)) +
  facet_wrap(~lab, scale = "free", ncol = 6) +
  ylab("Expression") + xlab(NULL) +
  ggtitle("GSE87466") +
  scale_color_bmj(name = "Group") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top", 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave("03.Model/07.GSE87466.Box.png", width = 8, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("03.Model/07.GSE87466.Box.pdf", width = 8, height = 4, units = "in", dpi = 300, bg = "white")



rm(list = ls())
df.exp = read.delim2("cleandata/GSE87466.exp.xls", row.names = 1) |> 
  type.convert(as.is = T)
df.pheno = read.delim2("cleandata/GSE87466.pheno.xls")
gs = read.delim2("02.Target/06.Intersection.xls")
df.exp = df.exp[gs$gene,df.pheno$geo_accession]

res.roc = apply(df.exp, 1, function(x){
  roc(df.pheno$group, x)
})

res.auc = sapply(res.roc, pROC::auc)
res.auc = data.frame(Model = names(res.auc), AUC = res.auc)
res.auc$AUC = round(res.auc$AUC, 3)
res.auc$lab = paste0(res.auc$Model, " (AUC = ", res.auc$AUC, ")")
res.auc = res.auc[order(res.auc$AUC, decreasing = T),]
res.roc = res.roc[res.auc$Model]
names(res.roc) = res.auc$lab

ggroc(res.roc, legacy.axes = T, size = 0.8) + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "grey50", linetype = 2) +
  theme_bw() + 
  ggtitle("GSE87466") +
  theme(panel.border = element_rect(color = "black", fill = "transparent"), 
        aspect.ratio = 1, legend.title = element_blank(), 
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"), 
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
ggsave("03.Model/08.GSE87466.ROC.png", width = 7, height = 5.5, units = "in", dpi = 300, bg = "white")
ggsave("03.Model/08.GSE87466.ROC.pdf", width = 7, height = 5.5, units = "in", dpi = 300, bg = "white")



rm(list = ls())
df.exp = read.delim2("cleandata/GSE59071.exp.xls", row.names = 1) |> 
  type.convert(as.is = T)
df.pheno = read.delim2("cleandata/GSE59071.pheno.xls")
gs = read.delim2("02.Target/06.Intersection.xls")
df.exp = df.exp[gs$gene,df.pheno$geo_accession]

p.diff = apply(df.exp, 1, \(x){
  wilcox.test(x~df.pheno$group)$p.value
})
df.diff = data.frame(gene = rownames(df.exp),
                     p.val = p.diff,
                     p.adj = p.adjust(p.diff, method = "BH"))
df.diff$lab = paste0(df.diff$gene, "\np.adj = ", signif(df.diff$p.adj,3))

df.exp$gene = rownames(df.exp)
df.exp = reshape2::melt(df.exp, id.var = "gene", variable.name = "sample")
df.exp$group = df.pheno$group[match(df.exp$sample, df.pheno$geo_accession)]
df.exp$lab = df.diff$lab[match(df.exp$gene, df.diff$gene)]

ggplot(df.exp, aes(x = group, y = value)) + 
  geom_boxplot(aes(color = group)) +
  facet_wrap(~lab, scale = "free", ncol = 6) +
  ylab("Expression") + xlab(NULL) +
  ggtitle("GSE59071") +
  scale_color_bmj(name = "Group") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top", 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave("03.Model/09.GSE59071.Box.png", width = 8, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("03.Model/09.GSE59071.Box.pdf", width = 8, height = 4, units = "in", dpi = 300, bg = "white")



rm(list = ls())
df.exp = read.delim2("cleandata/GSE59071.exp.xls", row.names = 1) |> 
  type.convert(as.is = T)
df.pheno = read.delim2("cleandata/GSE59071.pheno.xls")
gs = read.delim2("02.Target/06.Intersection.xls")
df.exp = df.exp[gs$gene,df.pheno$geo_accession]

res.roc = apply(df.exp, 1, function(x){
  roc(df.pheno$group, x)
})

res.auc = sapply(res.roc, pROC::auc)
res.auc = data.frame(Model = names(res.auc), AUC = res.auc)
res.auc$AUC = round(res.auc$AUC, 3)
res.auc$lab = paste0(res.auc$Model, " (AUC = ", res.auc$AUC, ")")
res.auc = res.auc[order(res.auc$AUC, decreasing = T),]
res.roc = res.roc[res.auc$Model]
names(res.roc) = res.auc$lab

ggroc(res.roc, legacy.axes = T, size = 0.8) + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "grey50", linetype = 2) +
  theme_bw() + 
  ggtitle("GSE59071") +
  theme(panel.border = element_rect(color = "black", fill = "transparent"), 
        aspect.ratio = 1, legend.title = element_blank(), 
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"), 
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
ggsave("03.Model/10.GSE59071.ROC.png", width = 7, height = 5.5, units = "in", dpi = 300, bg = "white")
ggsave("03.Model/10.GSE59071.ROC.pdf", width = 7, height = 5.5, units = "in", dpi = 300, bg = "white")




