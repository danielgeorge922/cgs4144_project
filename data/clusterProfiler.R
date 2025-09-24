library(readr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(rstudioapi)

data <- read_csv("significant_genes.csv")
genes <- data %>%
  filter(`p-value` < 0.05) %>%
  pull(GeneName)
genes <- as.character(genes)

entrez_genes <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
entrez_ids <- entrez_genes$ENTREZID

ego <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",        # Change to "MF" or "CC" if desired
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

head(ego)
write.csv(as.data.frame(ego), "GO_enrichment_results.csv", row.names = FALSE)
dotplot(ego) + ggtitle("GO Biological Process Enrichment")
