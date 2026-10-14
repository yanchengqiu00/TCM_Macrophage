rm(list = ls())
library(org.Hs.eg.db)
library(tidyverse)
library(GEOquery)

gs = getGEO(filename = "rawdata/GSE87466_series_matrix.txt.gz", getGPL = F)
df.exp = Biobase::exprs(gs) |> as.data.frame()
df.exp = limma::normalizeBetweenArrays(df.exp) |> as.data.frame()
df.pheno = gs |> pData() |>  as.data.frame()
df.pheno$sample = df.pheno$geo_accession
df.pheno$group = ifelse(df.pheno$`disease:ch1`=="Ulcerative Colitis (UC)", "Disease","Control")
df.annot = getGEO(filename = "rawdata/GPL13158.annot.gz") |> Table()
df.annot = subset(df.annot, !str_detect(`Gene symbol`, "//"))
df.exp$gene = df.annot$`Gene symbol`[match(rownames(df.exp), df.annot$ID)]
df.exp = subset(df.exp, gene != "")
df.exp = na.omit(df.exp)
df.exp = aggregate(.~gene, FUN = max, data = df.exp)
write.table(df.exp, "cleandata/GSE87466.exp.xls", sep = "\t", quote = F, row.names = F)
write.table(df.pheno, "cleandata/GSE87466.pheno.xls", sep = "\t", quote = F, row.names = F)


rm(list = ls())
gs = getGEO(filename = "rawdata/GSE59071_series_matrix.txt.gz", getGPL = F)
df.exp = gs |>  exprs() |>  as.data.frame()
df.pheno = gs |>  pData() |>  as.data.frame()
df.annot = getGEO(filename = "rawdata/GPL6244.annot.gz") |> Table()
df.exp$gene = df.annot$`Gene symbol`[match(rownames(df.exp), df.annot$ID)]
df.exp = subset(df.exp, gene != "")
df.exp = na.omit(df.exp)
df.exp = df.exp[!str_detect(df.exp$gene, fixed("/")),]
df.exp = aggregate(.~gene, FUN = max, data = df.exp)
df.pheno$sample = df.pheno$geo_accession
df.pheno = df.pheno[df.pheno$`disease:ch1` != "Crohn's disease",]
df.pheno = df.pheno[df.pheno$`disease activity:ch1` != "inactive",]
df.exp = df.exp[c("gene",df.pheno$sample)]
df.pheno$group = ifelse(str_detect(df.pheno$title, "control"), "Control", "Disease")
write.table(df.exp, "cleandata/GSE59071.exp.xls", sep = "\t", quote = F, row.names = F)
write.table(df.pheno, "cleandata/GSE59071.pheno.xls", sep = "\t", quote = F, row.names = F)


