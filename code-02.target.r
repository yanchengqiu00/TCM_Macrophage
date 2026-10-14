rm(list = ls())
library(rvest)
library(httr)
library(jsonlite)
library(dplyr)
library(data.table)
library(tidyverse)
library(org.Hs.eg.db)

obThresholdValue = 30
dlThresholdValue = 0.18

names <- c(
  "Aucklandiae Radix", #muxiang
  "Coptidis Rhizoma", #huanglian
  "Bletilla Striata", #baiji
  "Radix Sanguisorbae", #diyu
  "Agrimonia Eupatoria", #xianhecao
  "Cortex Moutan", #mudanpi
  "Licorice", #gancao
  "Scutellariae Radix", #huangqin
  "Radix Puerariae" #gegen
)

urls <- c(
  "https://www.tcmsp-e.com/tcmspsearch.php?qr=Aucklandiae%20Radix&qsr=herb_en_name&token=be5636c331b0a0e0c97dae71d5cfab0f",
  "https://www.tcmsp-e.com/tcmspsearch.php?qr=Coptidis%20Rhizoma&qsr=herb_en_name&token=be5636c331b0a0e0c97dae71d5cfab0f",
  "https://www.tcmsp-e.com/tcmspsearch.php?qr=Bletilla%20Striata%20(Thunb.Ex%20A.Murray)Rchb.F.&qsr=herb_en_name&token=be5636c331b0a0e0c97dae71d5cfab0f",
  "https://www.tcmsp-e.com/tcmspsearch.php?qr=Radix%20Sanguisorbae&qsr=herb_en_name&token=be5636c331b0a0e0c97dae71d5cfab0f",
  "https://www.tcmsp-e.com/tcmspsearch.php?qr=Agrimonia%20Eupatoria&qsr=herb_en_name&token=be5636c331b0a0e0c97dae71d5cfab0f",
  "https://www.tcmsp-e.com/tcmspsearch.php?qr=Cortex%20Moutan&qsr=herb_en_name&token=be5636c331b0a0e0c97dae71d5cfab0f",
  "https://www.tcmsp-e.com/tcmspsearch.php?qr=licorice&qsr=herb_en_name&token=be5636c331b0a0e0c97dae71d5cfab0f",
  "https://www.tcmsp-e.com/tcmspsearch.php?qr=Scutellariae%20Radix&qsr=herb_en_name&token=be5636c331b0a0e0c97dae71d5cfab0f",
  "https://www.tcmsp-e.com/tcmspsearch.php?qr=Radix%20Puerariae&qsr=herb_en_name&token=be5636c331b0a0e0c97dae71d5cfab0f"
)

uniprotData <-
  tibble(fread("../database/uniprotkb_AND_model_organism_9606_2024_11_11.tsv.gz")) %>%
  filter(!is.na(`Gene Names`)) %>%
  dplyr::select(`Protein names`,`Gene Names`) %>%
  mutate(`Gene Names` = gsub('\\s.*', '', `Gene Names`), Symbol = gsub('\\s.*', '', `Gene Names`)) %>%
  filter(!str_detect(Symbol, '\\[a-z\\]|\\-|\\_')) %>%
  mutate(
    protein_names_space = gsub('-', ' ', `Protein names`),
    protein_names_first = gsub('\\s+\\(.*', '', `Protein names`),
    protein_names_first_space = gsub('-', ' ', gsub('\\s+\\(.*', '', `Protein names`))
  )
symbol_all <- tibble(Drug=NA,MOL_ID=NA,molecule_name=NA,target_name=NA,Symbol=NA)

num <- 1
select = dplyr::select
for (url in urls) {
  name <- names[num]
  print(name)
  num <- num + 1
  web <- read_html(GET(url,encoding="UTF-8", config(ssl_verifypeer = FALSE)))
  
  tcmsp <- web %>% html_elements("script") %>% html_text()
  test1 <- str_extract_all(tcmsp,"data:\\s\\[.*\\]")
  test2 <- unlist(test1[12])
  
  drug_m <- str_replace(test2[1], "data:","") %>% fromJSON(simplifyVector = TRUE) %>% mutate(across(c(ob, dl), as.numeric)) %>% filter(ob >= obThresholdValue & dl >= dlThresholdValue)
  drug_t <- str_replace(test2[2], "data:","") %>% fromJSON(simplifyVector = TRUE) %>% semi_join(drug_m, by = "MOL_ID") %>% tibble() %>% add_column(Drug = name, .before = 'molecule_ID')
  #drug_t = subset(drug_t, validated = "validated")
  all_symbol <- drug_t %>% left_join(uniprotData, join_by(target_name == `Protein names`)) %>% filter(!is.na(`Gene Names`)) %>% select(Drug,MOL_ID, molecule_name, target_name, Symbol)
  all_symbol1 <- drug_t %>% left_join(uniprotData, join_by(target_name == protein_names_space)) %>% filter(!is.na(`Gene Names`)) %>% select(Drug,MOL_ID, molecule_name, target_name, Symbol)
  all_symbol2 <- drug_t %>%left_join(uniprotData, join_by(target_name == protein_names_first)) %>% filter(!is.na(`Gene Names`)) %>% select(Drug,MOL_ID, molecule_name, target_name, Symbol)
  all_symbol3 <- drug_t %>%left_join(uniprotData, join_by(target_name == protein_names_first_space)) %>% filter(!is.na(`Gene Names`)) %>% select(Drug,MOL_ID, molecule_name, target_name, Symbol)
  
  symbol_all <- symbol_all %>% add_row(all_symbol) %>% add_row(all_symbol1) %>% add_row(all_symbol2) %>% add_row(all_symbol3) %>% distinct(target_name, molecule_name,Drug, .keep_all = TRUE)
}
symbol_all <- symbol_all %>% filter(!is.na(Symbol)) %>% dplyr::rename(MolId=MOL_ID,MolName=molecule_name,TargetName=target_name)
symbol_all = subset(symbol_all, Symbol != "")

df.annot = AnnotationDbi::select(org.Hs.eg.db, unique(symbol_all$Symbol), "ENTREZID", "SYMBOL")
symbol_all$Symbol = ifelse(symbol_all$Symbol == "opn", "SPP1", symbol_all$Symbol)
symbol_all$Symbol = ifelse(symbol_all$Symbol == "carbonic", "CA2", symbol_all$Symbol)
symbol_all$Symbol = ifelse(symbol_all$Symbol == "TNFA", "TNF", symbol_all$Symbol)
symbol_all$Symbol = ifelse(symbol_all$Symbol == "ACC1", "ACACA", symbol_all$Symbol)

symbol_all = symbol_all[c(1,3,5)]
colnames(symbol_all) = c("Drug","Molecule","Gene")
write.table(symbol_all, "02.Target/01.Drug.Targets.xls", sep = "\t", row.names = F, quote = F)



rm(list = ls())
df.tar = read.delim2("02.Target/01.Drug.Targets.xls")
df.1 = df.tar[1:2]
df.2 = df.tar[2:3]
colnames(df.1) = colnames(df.2) = c("source","target")
df.net = rbind(df.1, df.2)
df.net = df.net[!duplicated(df.net),]
write.table(df.net, "interm/all.net.txt", quote = F, sep = "\t", row.names = F)
write.table(data.frame(unique(df.tar$Molecule)), "interm/mols.txt", quote = F, sep = "\t", row.names = F)
write.table(data.frame(unique(df.tar$Gene)), "interm/genes.txt", quote = F, sep = "\t", row.names = F)





rm(list = ls())
library(tidyverse)
library(limma)
library(ggrepel)

df = read.delim2("cleandata/GSE87466.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
df.pheno = read.delim2("cleandata/GSE87466.pheno.xls")
df = df[df.pheno$geo_accession]
design.mat = cbind(Control = ifelse(df.pheno$group == "Control", 1, 0), 
                   Patient = ifelse(df.pheno$group == "Control", 0, 1))
contrast.mat = makeContrasts(contrasts="Patient-Control", levels=design.mat)

fit = lmFit(df, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
fit = fit[c(1,4,5)]
fit = cbind(Gene = rownames(fit), fit)
fit$Direction = ifelse(fit$logFC > 0, "Up", "No")
fit$Direction = ifelse(fit$logFC < 0, "Down", fit$Direction)
fit$Direction = ifelse(fit$adj.P.Val < 0.05, fit$Direction, "No")
fit$Direction = factor(fit$Direction, levels = c("Down", "No", "Up"))
fit = fit[order(fit$adj.P.Val, decreasing = F), ]
fit.out = subset(fit, Direction != "No")
write.table(fit.out, "02.Target/03.DEG.xls", sep = "\t", quote = F, col.names = T, row.names = F)

top = rbind(subset(fit.out, Direction == "Up")[1:5,],
            subset(fit.out, Direction == "Down")[1:5,])

ggplot(fit, aes(x = logFC, y = -log10(adj.P.Val), color = Direction)) + 
  geom_point(size = 2, alpha = 0.7) + 
  geom_label_repel(aes(label = Gene), data = top, show.legend = F) +
  scale_color_manual(values=c("darkblue", "grey30", "darkred"), labels = c("Down", "No Change", "Up")) +
  geom_hline(yintercept=-log10(0.05), col="grey60", linetype = 2) +
  ggtitle("Disease vs Control") +
  ylab("-log10(p.adj)") +
  theme_minimal(base_size = 16) + 
  theme(legend.title = element_blank(), 
        legend.position = "top", 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5) )
ggsave("02.Target/04.Volcano.png", width = 10, height = 6, dpi = 300, bg = "white")
ggsave("02.Target/04.Volcano.pdf", width = 10, height = 6, dpi = 300, bg = "white")



rm(list = ls())
library(ggvenn)
df.deg = read.delim2("02.Target/03.DEG.xls")
df.deg.2 = read.delim2("../20250610-CCPJ-2002-5/02.Target/02.DEG.xls")
df.tar = read.delim2("02.Target/01.Drug.Targets.xls")
df.sc = read.delim2("01.SC/14.Gene.Module.xls")
df.sc = df.sc[df.sc$module %in% c("turquoise"),]

ggvenn(list(DEG = df.deg$Gene, Target = df.tar$Gene, SC = df.sc$gene_name), 
       stroke_size = 0.3, set_name_size = 6, text_size = 5, show_percentage = F) + 
  ylim(c(-2,2)) + 
  scale_fill_bmj()
ggsave("02.Target/05.Venn.png", width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("02.Target/05.Venn.pdf", width = 5, height = 5, units = "in", dpi = 300, bg = "white")

inter = intersect(df.deg$Gene, df.tar$Gene) |> intersect(df.sc$gene)
write.table(data.frame(gene = inter), "02.Target/06.Intersection.xls", quote = F, sep = "\t", row.names = F)

df.tar = df.tar[df.tar$Gene %in% inter,]
df.drug = df.tar[1:2]
df.gene = df.tar[2:3]
df.drug = df.drug[!duplicated(df.drug),]
df.gene = df.gene[!duplicated(df.gene),]
colnames(df.drug) = colnames(df.gene) = c("source","target")
df.out = rbind(df.drug, df.gene)
write.table(df.out, "interm/inter.mol.drug.txt", quote = F, sep = "\t", row.names = F)
write.table(data.frame(drug = unique(df.drug$source)), "interm/drugs.txt", quote = F, sep = "\t", row.names = F)
write.table(data.frame(mol = unique(df.drug$target)), "interm/mols.txt", quote = F, sep = "\t", row.names = F)
write.table(data.frame(gene = unique(df.gene$target)), "interm/genes.txt", quote = F, sep = "\t", row.names = F)



rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(gground)
df.deg = read.delim2("02.Target/06.Intersection.xls")
gs = unique(df.deg$gene)
df.annot = AnnotationDbi::select(org.Hs.eg.db, gs, "PATH", "SYMBOL")
df.annot$PATH = paste0("ko", df.annot$PATH)
df.kegg = read.delim2("../database/KEGG/kegg.txt")
df.annot$pathway = df.kegg$L3[match(df.annot$PATH, df.kegg$PathwayID)]
df.annot = df.annot[c(1,3)] |> na.omit()

df.tar = read.delim2("02.Target/01.Drug.Targets.xls")
df.tar = df.tar[df.tar$Gene %in% gs,]
df.1 = df.tar[1:2]
df.2 = df.tar[2:3]
colnames(df.1) = colnames(df.2) = colnames(df.annot) = c("source","target")
df.net = rbind(df.1, df.2, df.annot)
df.net = df.net[!duplicated(df.net),]
write.table(df.net, "interm/hub.net.txt", quote = F, sep = "\t", row.names = F)


