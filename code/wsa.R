library(Seurat)
setwd("/cluster/home/ztao_jh/projects/anteater/data/schuman/yixianai")

# objfiles <- list.files("./merged/merge100/", full.names = T)
# objlist <- parallel::mclapply(objfiles, readRDS, mc.cores = 13L)
# objlist <- parallel::mclapply(objlist, function(x){
#   cn <- str_remove(colnames(x@meta.data), pattern = "id_[0-9_]+")
#   colnames(x@meta.data) <- cn
#   x
# }, mc.cores = 13L)
# mergeobj <- merge(objlist[[1]], objlist[-1])
# 
# mergeobj <- NormalizeData(mergeobj, normalization.method = "CLR") %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# v <- c("WSI_normal_mask_full_pix", "close_margin_mask_full_pix", "WSI_tumor_mask_full_pix")
# lb <- apply(mergeobj@meta.data[v], 1, function(x) { v[which.max(x)] })
# mergeobj@meta.data$lable <- lb
# mergeobj@meta.data$lable_full <- paste0(str_remove(mergeobj@meta.data$orig.ident, pattern = "pix"),"_",lb)
# library(harmony)
# mergeobj <- RunHarmony(mergeobj, "orig.ident")
# mergeobj <- RunUMAP(mergeobj,  dims = 1:30, reduction = "harmony")
# pdf("umap100_harmony.pdf")
# DimPlot(mergeobj, reduction = "umap", label=T, group.by = "orig.ident")
# FeaturePlot(mergeobj, reduction = "umap", features = "WSI_normal_mask_full_pix", order = T)
# FeaturePlot(mergeobj, reduction = "umap", features = "close_margin_mask_full_pix", order = T)
# FeaturePlot(mergeobj, reduction = "umap", features = "WSI_tumor_mask_full_pix", order = T)
# DimPlot(mergeobj, reduction = "umap", group.by = "lable", order = F) + theme(legend.position = "top")
# DimPlot(mergeobj, reduction = "umap", group.by = "lable_full", order = F) + 
#   theme(legend.position = "top") + scale_color_manual(
#     values = c(
#       "#FFA07A","#FA8072","#F08080",
#                "#E9967A","#CD5C5C","#8B0000","#FFD700","#FFA500",
#                "#FF8C00","#EEE8AA","#F0E68C","#BDB76B","#7FFF00","#32CD32","#006400","#808000",
#                "#556B2F","#6B8E23","#00FF7F","#00FA9A","#8FBC8F","#E0FFFF","#00FFFF","#7FFFD4",
#                "#5F9EA0","#40E0D0","#008B8B","#00CED1","#20B2AA","#008080","#B0E0E6","#ADD8E6",
#                "#87CEFA","#D8BFD8","#DDA0DD","#DA70D6","#FFB6C1","#FF1493","#C71585"
#     )
#   )
# dev.off()
# 
# write_rds(mergeobj, file = "./merged/merge100.rds")
mergeobj <- "./merged/merge100.rds" %>% read_rds()
pdf("umap100_harmony0.pdf")
DimPlot(mergeobj, reduction = "umap", label=T, group.by = "orig.ident") + NoLegend()
DimPlot(mergeobj, reduction = "umap", group.by = "lable", order = F) + theme(legend.position = "top") + NoLegend()
DimPlot(mergeobj, reduction = "umap", group.by = "lable_full", order = F) + 
  theme(legend.position = "top") + scale_color_manual(
    values = c(
      "#FFA07A","#FA8072","#F08080",
               "#E9967A","#CD5C5C","#8B0000","#FFD700","#FFA500",
               "#FF8C00","#EEE8AA","#F0E68C","#BDB76B","#7FFF00","#32CD32","#006400","#808000",
               "#556B2F","#6B8E23","#00FF7F","#00FA9A","#8FBC8F","#E0FFFF","#00FFFF","#7FFFD4",
               "#5F9EA0","#40E0D0","#008B8B","#00CED1","#20B2AA","#008080","#B0E0E6","#ADD8E6",
               "#87CEFA","#D8BFD8","#DDA0DD","#DA70D6","#FFB6C1","#FF1493","#C71585"
    )
  ) + NoLegend()
dev.off()
# ==============================================================================
library(Seurat)
options(digits = 15)
setwd("/cluster/home/ztao_jh/projects/anteater/data/schuman/yixianai")
# objfiles <- list.files("./merged/merge200/", full.names = T, pattern = ".rds")
# objlist <- parallel::mclapply(objfiles, readRDS, mc.cores = 13L)
# 
# df1 <- "./merged/merge200/80551-2.csv" %>% read.csv() %>% # 2
#   dplyr::filter(error < 0.1)
# df2 <- "./merged/merge200/82020.csv" %>% read.csv() %>%  # 11
#   dplyr::filter(error < 0.1)
# 
# cof <- base::intersect(df1$ionx0, df2$ionx0)
# df1 <- df1 %>% dplyr::filter(ionx0 %in% {{ cof }})
# df2 <- df2 %>% dplyr::filter(ionx0 %in% {{ cof }})
# 
# mx2 <- objlist[[2]]@assays$MS@counts[df1$ionx,]
# rownames(mx2) <- df1$ionx0
# df <- objlist[[2]]@meta.data
# cn <- str_remove(colnames(df), pattern = "id_[0-9_]+")
# colnames(df) <- cn
# obj2 <- CreateSeuratObject(counts = mx2, meta.data = df)
# 
# mx12 <- objlist[[11]]@assays$MS@counts[df2$ionx,]
# df <- objlist[[11]]@meta.data
# cn <- str_remove(colnames(df), pattern = "id_[0-9_]+")
# colnames(df) <- cn
# rownames(mx12) <- df2$ionx0
# obj12 <- CreateSeuratObject(counts = mx12, meta.data = df)
# 
# objlist0 <- parallel::mclapply(objlist[-c(2,11)], function(x){
#   df <- x@meta.data
#   cn <- str_remove(colnames(df), pattern = "id_[0-9_]+")
#   colnames(df) <- cn
# 
#   CreateSeuratObject(counts = x@assays$MS@counts, meta.data = df)
# }, mc.cores = 13L)
# objlistx <- c(objlist0, obj2, obj12)
# 
# 
# mergeobj <- merge(objlistx[[1]], objlistx[-1])
# zn <- apply(mergeobj@assays$RNA@counts, 1, function(x) sum(x == 0))
# zn %>% density() %>% plot()
# mergeobj <- mergeobj[zn <= 3000,]
# mergeobj <- NormalizeData(mergeobj, normalization.method = "CLR") %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# 
# v <- c("WSI_normal_mask_full_pix", "close_margin_mask_full_pix", "WSI_tumor_mask_full_pix")
# lb <- apply(mergeobj@meta.data[v], 1, function(x) { v[which.max(x)] })
# mergeobj@meta.data$lable <- lb
# mergeobj@meta.data$lable_full <- paste0(str_remove(mergeobj@meta.data$orig.ident, pattern = "pix"),"_",lb)
# library(harmony)
# mergeobj <- RunHarmony(mergeobj, "orig.ident")
# mergeobj <- RunUMAP(mergeobj,  dims = 1:30, reduction = "harmony")
# 
# pdf("umap200_harmony.pdf")
# DimPlot(mergeobj, reduction = "umap", label=T, group.by = "orig.ident")
# FeaturePlot(mergeobj, reduction = "umap", features = "WSI_normal_mask_full_pix", order = T)
# FeaturePlot(mergeobj, reduction = "umap", features = "close_margin_mask_full_pix", order = T)
# FeaturePlot(mergeobj, reduction = "umap", features = "WSI_tumor_mask_full_pix", order = T)
# DimPlot(mergeobj, reduction = "umap", group.by = "lable", order = F) + theme(legend.position = "top")
# DimPlot(mergeobj, reduction = "umap", group.by = "lable_full", order = F) + 
#   theme(legend.position = "top") + 
#   scale_color_manual(
#   values = c(
#     "#FFA07A","#FA8072","#F08080",
#              "#E9967A","#CD5C5C","#8B0000","#FFD700","#FFA500",
#              "#FF8C00","#EEE8AA","#F0E68C","#BDB76B","#7FFF00","#32CD32","#006400","#808000",
#              "#556B2F","#6B8E23","#00FF7F","#00FA9A","#8FBC8F","#E0FFFF","#00FFFF","#7FFFD4",
#              "#5F9EA0","#40E0D0","#008B8B","#00CED1","#20B2AA","#008080","#B0E0E6","#ADD8E6",
#              "#87CEFA","#D8BFD8","#DDA0DD","#DA70D6","#FFB6C1","#FF1493","#C71585"
#   )
# )
# dev.off()
# ncol(mergeobj)
# write_rds(mergeobj, file = "./merged/merge200.rds")


mergeobj <- "./merged/merge200.rds" %>% read_rds()
pdf("umap200_harmony0.pdf")
DimPlot(mergeobj, reduction = "umap", label=T, group.by = "orig.ident") + NoLegend()
DimPlot(mergeobj, reduction = "umap", group.by = "lable", order = F) + theme(legend.position = "top") + NoLegend()
DimPlot(mergeobj, reduction = "umap", group.by = "lable_full", order = F) + 
  theme(legend.position = "top") + scale_color_manual(
    values = c(
      "#FFA07A","#FA8072","#F08080",
               "#E9967A","#CD5C5C","#8B0000","#FFD700","#FFA500",
               "#FF8C00","#EEE8AA","#F0E68C","#BDB76B","#7FFF00","#32CD32","#006400","#808000",
               "#556B2F","#6B8E23","#00FF7F","#00FA9A","#8FBC8F","#E0FFFF","#00FFFF","#7FFFD4",
               "#5F9EA0","#40E0D0","#008B8B","#00CED1","#20B2AA","#008080","#B0E0E6","#ADD8E6",
               "#87CEFA","#D8BFD8","#DDA0DD","#DA70D6","#FFB6C1","#FF1493","#C71585"
    )
  ) + NoLegend()
dev.off()







