# based on https://alexslemonade.github.io/refinebio-examples/03-rnaseq/gene-id-annotation_rnaseq_01_ensembl.html
if (!("org.Hs.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Hs.eg.db", update = FALSE) # Use package for humans
}

library(org.Hs.eg.db)
library(magrittr)
library(readr)
library(dplyr)
library(tibble)

# Load data
metadata <- read_tsv("data/SRP075806/metadata_SRP075806.tsv")

expression_df <- read_tsv("data/SRP075806/SRP075806.tsv") %>%
  column_to_rownames("Gene") %>%
  select(all_of(metadata$refinebio_accession_code))

# Check ordering
stopifnot(all.equal(colnames(expression_df), metadata$refinebio_accession_code))

# Map Ensembl IDs to gene symbols
# Strip version numbers first (e.g. ENSG00000141510.15 -> ENSG00000141510)
ensembl_ids <- sub("\\..*$", "", rownames(expression_df))

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys     = ensembl_ids,
  keytype  = "ENSEMBL",
  column   = "SYMBOL",
  multiVals = "first"
)

# Add gene symbols & save
expression_df$Gene <- gene_symbols
expression_df <- expression_df %>% relocate(GeneSymbol, .before = 1)

# Drop rows without a mapped symbol if desired
expression_clean <- expression_df[!is.na(expression_df$GeneSymbol), ]

write_tsv(
  expression_clean,
  "data/SRP075806/SRP075806_gene_symbols.tsv"
)
