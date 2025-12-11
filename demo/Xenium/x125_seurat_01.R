
output_dir <- "/home/zhangfeng/projects/QRLu/data/12.5/"

setwd(output_dir)
library(Seurat)
library(dplyr)
library(data.table)

count_matrix_dt <- fread(paste0(output_dir,'/',"xenium_count_matrix.csv"))
count_matrix <- as.matrix(count_matrix_dt, rownames = "cell_id")
count_matrix <- t(count_matrix)
cell_metadata <- read.csv(paste0(output_dir,'/',"xenium_cell_metadata.csv"),row.names = 1)

saveRDS(cell_metadata,'cell_metadata.rds')
saveRDS(count_matrix, 'count_matrix.rds')


library(Seurat)
library(dplyr)
counts = count_matrix[,2:ncol(count_matrix)]
metadata =cell_metadata[2:nrow(cell_metadata),]

seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
seurat_obj[["spatial"]] = CreateAssayObject(counts = counts)
DefaultAssay(seurat_obj) <- "spatial"

seurat_obj@images <- list()
coords=metadata[, c("x_centroid", "y_centroid")]
#colnames(coords)=c('centroid_x','centroid_y')

seurat_obj@images[["image"]] <- new(
  Class = "SlideSeq",  # 若为Visium数据则改为"Visium"
  assay = "spatial",   # 与DefaultAssay一致
  coordinates = coords,
  key = "image_"       # 用于关联元数据的前缀（必须）
)


pdf('SOX11.pdf')
SpatialFeaturePlot(seurat_obj, features = "SOX11")
dev.off()


saveRDS(seurat_obj,'seurat_obj_raw.rds')

###########

















