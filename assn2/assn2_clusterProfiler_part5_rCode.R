library(readr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(stringr)
library(tibble)

data <- read_csv("significant_genes.csv")

genes <- as.character(data$GeneName)

entrez_genes <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
entrez_ids <- entrez_genes$ENTREZID

write_csv(entrez_genes, "mapped_entrez_ids.csv")

ego <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

ego_df <- as.data.frame(ego)
write_csv(ego_df, "GO_enrichment_results_clusterprofiler.csv")

ego_table <- as_tibble(ego_df) %>%
  mutate(
    Gene_set = "GO_Biological_Process_2023",
    Term = paste0(Description, " (", ID, ")"),
    Overlap = paste0(Count, "/", as.character(str_extract(GeneRatio, "\\d+$"))),
    `P-value` = pvalue,
    `Adjusted P-value` = p.adjust,
    `Old P-value` = pvalue,
    `Old Adjusted P-value` = p.adjust,
    `Odds Ratio` = NA,
    `Combined Score` = NA,
    Genes = gsub("/", ";", geneID)
  )

ego_table <- dplyr::select(
  ego_table,
  Gene_set, Term, Overlap, `P-value`, `Adjusted P-value`,
  `Old P-value`, `Old Adjusted P-value`, `Odds Ratio`, `Combined Score`, Genes
)

write_csv(ego_table, "GO_enrichment_results_full_table.csv")

p <- barplot(ego, showCategory = 27, horiz = TRUE) +
  ggtitle("GO Biological Process Enrichment") +
  theme(axis.text.y = element_text(size = 8))
print(p)
