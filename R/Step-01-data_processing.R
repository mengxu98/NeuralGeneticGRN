library(dplyr)
library(tidyr)

load("../../data/stab_combined/stab_combined.Rdata")
astro4MetaData <- h.combined@meta.data %>% filter(cluster1 == "Astro1")
astro4 <- h.combined[, rownames(astro4MetaData)]
astro4Matrix <- astro4@assays$RNA@counts

save(astro4Matrix, astro4MetaData, file = "astro4Matrix.rda")

science_HAR <- read.csv("science_HAR.csv")

data_cleaned <- science_HAR %>%
  filter(associated_genes != "NONE") %>%
  mutate(associated_genes = gsub(" \\([^\\)]+\\)", "", associated_genes))

harGeneList <- strsplit(as.vector(data_cleaned$associated_genes), ", ")
harGenes <- purrr::list_c(harGeneList)
astro4Matrix <- as.matrix(astro4Matrix) %>% t()

astro4Net <- inferCSN::inferCSN(astro4Matrix,
                                regulators = harGenes,
                                # targets = harGenes,
                                cores = 6,
                                verbose = TRUE)

save(astro4Net, file = "astro4Net.rda")


network.heatmap(astro4Net[1:1000, ], heatmapSize = 30)
