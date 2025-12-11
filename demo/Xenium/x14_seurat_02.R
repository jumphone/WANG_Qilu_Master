
output_dir <- "/home/zhangfeng/projects/QRLu/data/Xenium_14/"
setwd(output_dir)

library(Seurat)
library(dplyr)
library(future)
getOption("future.globals.maxSize")
options(future.globals.maxSize = 100 * 1024^3)


seurat_obj=readRDS('seurat_obj_raw.rds')

seurat_obj <- SCTransform(seurat_obj, assay = "spatial", verbose = TRUE)

pdf('SOX11.pdf')
SpatialFeaturePlot(seurat_obj, features = "SOX11")
dev.off()

seurat_obj <- RunPCA(seurat_obj, assay = "SCT", verbose = TRUE)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:30)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)

saveRDS(seurat_obj,'seurat_obj_umap.rds')


pdf('output.pdf')
DimPlot(seurat_obj, reduction = "umap")
SpatialDimPlot(seurat_obj, cols = "polychrome")
dev.off()
















