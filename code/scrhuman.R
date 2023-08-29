library(Seurat)
library(tidyverse)
rna_merged <- "~/projects/anteater/data/schuman/yixianai/merged/rna_merged.rds" %>% read_rds()
DefaultAssay(rna_merged) <- "RNA"
rna_merged <- NormalizeData(rna_merged) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
rna_merged@assays$RNA@data %>% max


