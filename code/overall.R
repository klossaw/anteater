library(tidyverse)
library(glue)
library(Seurat)
library(CellChat)
setwd("~/projects/anteater/data/schuman/yixianai/")
# ======================== 1 =====================
objmege <- "./merged/rna_bigcelltype_merged.rds" %>% read_rds()
objmege <- RunTSNE(objmege,  dims = 1:30, reduction = "harmony")
Seurat::DefaultAssay(objmege)

objmege$lable_big0 <- objmege$lable_big
objmege$lable_big0[objmege$lable_detail %in% c("DC1DC2", "DC_CLEC9A")] <- "DC"
objmege$lable_big0[objmege$lable_detail %in%
                              c("Mac_CX3CR1_F4/80+", "Mac_APOE", "Mac_FCN1",
                                "Mac_CPA1","Mac_MKI67","Mono_SELENOP")] <- "Mac"
objmege$lable_big0[objmege$lable_detail %in%
                              c("Mono_KLF6")] <- "Mono"

pdf("total_umap_cd169_cd209.pdf")
DimPlot(objmege, group.by = "lable_big0", label = T)
FeaturePlot(objmege, features = "SIGLEC1", order = T)
FeaturePlot(objmege, features = "CD209", order = T)
dev.off()
pdf("total_tsne_cd169_cd209.pdf")
DimPlot(objmege, group.by = "lable_big0", reduction = "tsne", label = T)
FeaturePlot(objmege, features = "SIGLEC1", reduction = "tsne", order = T)
FeaturePlot(objmege, features = "CD209", reduction = "tsne", order = T)
dev.off()


# ======================== 2 =====================
myeloid <- "./myeloid/total_myeloid.rds" %>% read_rds()
DefaultAssay(myeloid) <- "RNA"
myeloid$lable_big0 <- myeloid$lable_big
myeloid$lable_big0[myeloid$lable_big %in% c("DC1DC2", "DC_CLEC9A")] <- "DC"
myeloid$lable_big0[myeloid$lable_big %in%
                     c("Mac_CX3CR1_F4/80+", "Mac_APOE", "Mac_FCN1",
                       "Mac_CPA1","Mac_MKI67","Mono_SELENOP")] <- "Mac"
myeloid$lable_big0[myeloid$lable_big %in%
                     c("Mono_KLF6")] <- "Mono"


myeloid$lable_big
DotPlot(myeloid, group.by = "lable_big0",
        features = c(
          "CD163","CD68",
          "CD14", "FCGR3A","FCGR1A",
          "CD1A", "CD1B",
          "SPP1", "APOE",
          "SIGLEC1")) +
  coord_flip()
VlnPlot(myeloid, group.by = "lable_big0",
        features = c(
          "CD163","CD68",
          "CD14", "FCGR3A","FCGR1A",
          "CD1A", "CD1B",
          "SPP1", "APOE",
          "SIGLEC1"))

DotPlot(myeloid, group.by = "lable_big0",
        features = c(
          "CD163","CD68","ADGRE1",
          "CD14", "FCGR3A","FCGR1A",
          "CD1A", "CD1B","ITGAM",
          "SPP1", "APOE",
          "SIGLEC1", "CD209")) +
  coord_flip()

DimPlot(myeloid, group.by = "lable_big0", label = T)
FeaturePlot(myeloid, features = c("CD68","SIGLEC1","CD209"), order = T)

library(ggsignif)
pdf("./cd169cd209plus_expression.pdf")
data.frame(
  CD169 = myeloid@assays$RNA@data["SIGLEC1",],
  group = myeloid$lable,
  celltype = myeloid$lable_big0
) %>% dplyr::filter(celltype == "Mac") %>% dplyr::filter(CD169 > 0) %>%
  ggplot(aes(x = group, y = CD169, fill = group)) +
  geom_boxplot() +
  theme_classic() +
  geom_signif(
    comparisons = list(c("Ctl","Her"),c("Her","Idio"), c("Ctl", "Idio")),
    map_signif_level = FALSE,
    step_increase = c(0.2,0.2,0.2))
data.frame(
  CD209 = myeloid@assays$RNA@data["CD209",],
  group = myeloid$lable,
  celltype = myeloid$lable_big0
) %>% dplyr::filter(celltype == "Mac") %>% dplyr::filter(CD209 > 0) %>%
  ggplot(aes(x = group, y = CD209, fill = group)) +
  geom_boxplot() +
  theme_classic() +
  geom_signif(
    comparisons = list(c("Ctl","Her"),c("Her","Idio"), c("Ctl", "Idio")),
    map_signif_level = FALSE,
    step_increase = c(0.2,0.2,0.2))
data.frame(
  CD169 = myeloid@assays$RNA@data["SIGLEC1",],
  group = myeloid$lable,
  celltype = myeloid$lable_big0
) %>% dplyr::filter(celltype == "Mac") %>% dplyr::filter(CD169 > 0) %>%
  ggplot(aes(x = group, y = CD169, fill = group)) +
  geom_boxplot() +
  theme_classic() +
  geom_signif(
    comparisons = list(c("Ctl","Her"),c("Her","Idio"), c("Ctl", "Idio")),
    map_signif_level = T,
    step_increase = c(0.2,0.2,0.2))
data.frame(
  CD209 = myeloid@assays$RNA@data["CD209",],
  group = myeloid$lable,
  celltype = myeloid$lable_big0
) %>% dplyr::filter(celltype == "Mac") %>% dplyr::filter(CD209 > 0) %>%
  ggplot(aes(x = group, y = CD209, fill = group)) +
  geom_boxplot() +
  theme_classic() +
  geom_signif(
    comparisons = list(c("Ctl","Her"),c("Her","Idio"), c("Ctl", "Idio")),
    map_signif_level = T,
    step_increase = c(0.2,0.2,0.2))



data.frame(
  CD169 = myeloid@assays$RNA@data["SIGLEC1",],
  group = myeloid$lable,
  celltype = myeloid$lable_big
) %>% dplyr::filter(celltype == "Mac_APOE") %>%
  ggplot(aes(x = group, y = CD169, fill = group)) +
  geom_violin(scale = "width") +
  theme_classic()

data.frame(
  CD169 = myeloid@assays$RNA@data["SIGLEC1",],
  group = myeloid$lable,
  celltype = myeloid$lable_big
) %>% dplyr::filter(celltype == "Mac_APOE") %>%
  ggplot(aes(x = group, y = CD169, fill = group)) +
  geom_boxplot() +
  theme_classic()

table(myeloid$lable, myeloid@assays$RNA@data["SPP1",] > 0)

data.frame(
  group = myeloid$lable,
  celltype = myeloid$lable_big
) %>% table()

data.frame(
  CD209 = myeloid@assays$RNA@data["CD209",],
  group = myeloid$lable,
  celltype = myeloid$lable_big0
) %>% dplyr::filter(celltype == "Mac") %>%
  ggplot(aes(x = group, y = CD209, fill = group)) +
  geom_violin(scale = "width") +
  theme_classic()


data.frame(
  lable = myeloid$lable,
  SPP1 = factor(myeloid@assays$RNA@data["SPP1",] > 0, levels = c("TRUE", "FALSE"))) %>%
  as.data.frame() %>%
  ggplot(aes(x = lable, fill = SPP1)) +
  geom_bar(position = "fill",color = "black") +
  theme_classic()
