# Setting working directory:
setwd("..")
# Loading packages:
library(tidyverse)
library(aracne.networks)
select <- dplyr::select

# Loading count data:
counts_data <- read_csv("data/johnstone2018/GSE101745_genewise-count.csv") %>%
  separate_rows(Name, sep = "/")
# Converting ENSEMBL ID to SYMBOL:
# Loading ENSEMBL aliases:
ENSEMBL <- read_tsv("data/ENSEMBL/mart_export.txt", col_names = c("Name", "SYMBOL"), skip = 1)
# Joining with counts data:
ENSEMBL %>% 
  filter(Name %in% counts_data$Name) %>% 
  unique() -> counts_genes
# Counting matches by SYMBOL:
counts_genes %>% 
  select(SYMBOL) %>% 
  drop_na() %>% 
  unique() %>% 
  nrow()
### Recovered 38372 matches between ENSEMBL ID and SYMBOL

# Writing regulon (with aracne.networks):
# write.regulon(regulonbrca, "tables/regulonbrca.tsv")

# Loading regulon:
regulon <- read_tsv("tables/regulonbrca.tsv", col_types = "ccdd")
# Extracting all gene IDs:
regulon_genes <- c(regulon$Regulator, regulon$Target) %>%
  unique() %>%
  enframe(name = NULL, value = "GENEID")
# Converting ENTREZ ID to SYMBOL:
# Loading NCBI aliases:
NCBI <- read_tsv("data/NCBI/Homo_sapiens.gene_info", col_types = cols(.default = "c"))
str(NCBI)
### Type of gene. Keep only protein-coding?
NCBI$type_of_gene <- factor(NCBI$type_of_gene)
levels(NCBI$type_of_gene)
# Adding SYMBOL to regulon:
NCBI %>%
  select(GENEID = GeneID, SYMBOL = Symbol, type_of_gene) %>%
  left_join(regulon_genes, .) -> regulon_genes
summary(regulon_genes)
### Almost all remaining genes are protein-coding. Keep others?


# Assessing SYMBOL overlap between counts data and regulon:
length(intersect(counts_genes$SYMBOL, regulon_genes$SYMBOL))
### 19209
# Filtering by counts data:
actual_regulon <- filter(regulon_genes, SYMBOL %in% counts_genes$SYMBOL)
summary(actual_regulon)
### We lose 304 genes (73 protein-coding)



# Result
gene_ids <- counts_genes %>%
  inner_join(regulon_genes) %>%
  drop_na(SYMBOL)

write_csv(gene_ids, "tables/gene_ids.csv")
