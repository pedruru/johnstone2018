setwd("..")

library(tidyverse)
library(DESeq2)
library(viper)
source("scripts/tinkering_viper.R")

counts_data <- read_csv("data/johnstone2018/GSE101745_genewise-count.csv") 
gene_ids <- read_csv("tables/gene_ids.csv", col_types = cols(.default = "c"))

counts_matrix <- counts_data %>% 
  filter(!(str_detect(Name, "/"))) %>% 
  inner_join(gene_ids) %>% 
  select(-c(Name, GENEID, type_of_gene)) %>% 
  column_to_rownames(var = "SYMBOL") %>% 
  data.matrix()

metadata <- data.frame(line = factor(c(rep("NIC", 2), rep(c("LNA", "LM2", "HM"), each = 6))),
                       batch = factor(c(rep("new", 2), rep(c(rep("new", 3), rep("old", 3)), 3)))) %>% 
  `rownames<-`(., colnames(counts_matrix))

res <- DESeqDataSetFromMatrix(countData = counts_matrix[, 9:20],
                              colData = metadata[9:20, ],
                              design = ~ batch + line)
res$line <- relevel(res$line, ref = "LM2")
res <- DESeq(res)
res <- data.frame(results(res, name = "line_HM_vs_LM2")[, c("stat", "pvalue")])

# nullmodel_reduced <- DESeq2Null(x = counts_matrix[, 15:20],
#                                 y = counts_matrix[, 9:14],
#                                 colData = metadata[9:20, ],
#                                 design = ~ batch + line,
#                                 refcol = "line",
#                                 reflvl = "LM2",
#                                 resname = "line_HM_vs_LM2",
#                                 cores = 4)
# 
# write_rds(nullmodel_reduced, file = "tables/nullmodel_reduced.rds")

# nullmodel_full <- DESeq2Null(x = counts_matrix,
#                              y = counts_matrix,
#                              colData = metadata,
#                              design = ~ batch + line,
#                              refcol = "line",
#                              reflvl = "LM2",
#                              resname = "line_HM_vs_LM2",
#                              cores = 4)
# 
# write_rds(nullmodel_full, file = "tables/nullmodel_full.rds")

#### PLOTTING ####

setwd("..")

library(tidyverse)

nullmodel_reduced <- read_rds("tables/nullmodel_reduced.rds")

# Sample 20 genes
ranres <- res %>% 
  filter(!is.na(.)) %>% 
  .[sample(nrow(.), 20), ]

rannullpvalue <- nullmodel_reduced[["pvalue"]] %>% 
  .[rownames(ranres), ]

for (i in 1:nrow(ranres)) {
  png(paste0("plots/ranres/nullmodel_p-values_for_", row.names(ranres[i, ]), ".png"))
  print(
    ggplot(data.frame("pvalue" = rannullpvalue[i, ]), aes(x = pvalue)) +
      geom_histogram(bins = 100, color = "black", fill = "grey") +
      geom_vline(aes(xintercept = ranres[i, ]$pvalue, color = "vline")) +
      scale_color_manual(name = "", values = c(vline = "red"), label = c(vline = "DESeq2 result")) +
      xlab(paste0("nullmodel p-values for ", rownames(ranres[i, ]))) +
      theme_classic() +
      theme(legend.position = c(0.85, 0.95))
  )
  dev.off()
}

# Permutation test p-values
res_stat <- filter(res, !is.na(res$stat))
nullmodel_reduced_stat <- nullmodel_reduced[["stat"]][rownames(res_stat), ]

pvalues <- data.frame()
for (i in 1:nrow(res_stat)) {
  pvalues[i, 1] <- (length(which(abs(nullmodel_reduced_stat[i, ]) >= abs(res_stat[i, 1]))) + 1) / 
    (length(nullmodel_reduced_stat[i, ]) + 1)
}
rownames(pvalues) <- rownames(res_stat)

pvalues %>% 
  ggplot(aes(x = V1)) +
  geom_histogram(bins = 100, color = "black", fill = "grey") +
  xlab("permutation test p-values") +
  theme_classic()

cbind(pvalues, res_stat$pvalue) %>% 
  ggplot(aes(x = V1, y = res_stat$pvalue)) +
  geom_point(alpha = 0.1) +
  xlab("permutation test p-values") +
  ylab("DESeq2 result p-values") +
  theme_classic()

cbind(pvalues, res_stat$pvalue) %>% 
  ggplot(aes(x = -log10(V1), y = -log10(res_stat$pvalue))) +
  geom_point(alpha = 0.1) +
  xlab("-log10 permutation test p-values") +
  ylab("-log10 DESeq2 result p-values") +
  theme_classic()
