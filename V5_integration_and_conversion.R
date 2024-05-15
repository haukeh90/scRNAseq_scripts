library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(reticulate)
library(dplyr)
library(data.table)
use_condaenv('/Users/hauke/anaconda3/')

#Run One-line integration-------------------------------------------------------
seurat_object <- readRDS("~/directory/seurat_object.rds")
DefaultAssay(seurat_object) <- "RNA5"

seurat_object <- IntegrateLayers(
  object = seurat_object, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = T
)

seurat_object <- IntegrateLayers(
  object = seurat_object, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = T
)

seurat_object <- IntegrateLayers(
  object = seurat_object, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = T
)

seurat_object <- IntegrateLayers(
  object = seurat_object, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = T
)

seurat_object <- IntegrateLayers(
  object = seurat_object, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/Users/hauke/anaconda3/", verbose = T
)

seurat_object <- FindNeighbors(seurat_object, reduction = "integrated.cca", dims = 1:30)
seurat_object <- FindClusters(seurat_object, resolution = 0.5, cluster.name = "cca_clusters", algorithm = 4)
seurat_object <- FindNeighbors(seurat_object, reduction = "integrated.rpca", dims = 1:30)
seurat_object <- FindClusters(seurat_object, resolution = 0.5, cluster.name = "rpca_clusters", algorithm = 4)
seurat_object <- FindNeighbors(seurat_object, reduction = "harmony", dims = 1:30)
seurat_object <- FindClusters(seurat_object, resolution = 0.5, cluster.name = "harmony_clusters", algorithm = 4)
seurat_object <- FindNeighbors(seurat_object, reduction = "integrated.mnn", dims = 1:30)
seurat_object <- FindClusters(seurat_object, resolution = 0.5, cluster.name = "mnn_clusters", algorithm = 4)
seurat_object <- FindNeighbors(seurat_object, reduction = "integrated.scvi", dims = 1:30)
seurat_object <- FindClusters(seurat_object, resolution = 0.5, cluster.name = "scvi_clusters", algorithm = 4)

seurat_object <- RunUMAP(seurat_object, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
seurat_object <- RunUMAP(seurat_object, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
seurat_object <- RunUMAP(seurat_object, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
seurat_object <- RunUMAP(seurat_object, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")
seurat_object <- RunUMAP(seurat_object, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")

p1 <- DimPlot(
  seurat_object,
  reduction = "umap.cca",
  group.by = c("conditions", "cca_clusters"),
  combine = FALSE, label.size = 2, shuffle = T
) 
p2 <- DimPlot(
  seurat_object,
  reduction = "umap.rpca",
  group.by = c("conditions", "rpca_clusters"),
  combine = FALSE, label.size = 2
) 
p3 <- DimPlot(
  seurat_object,
  reduction = "umap.harmony",
  group.by = c("conditions", "harmony_clusters"),
  combine = FALSE, label.size = 2, shuffle = T
) 
p4 <- DimPlot(
  seurat_object,
  reduction = "umap.mnn",
  group.by = c("conditions", "mnn_clusters"),
  combine = FALSE, label.size = 2, shuffle = T
) 

p5 <- DimPlot(
  seurat_object,
  reduction = "umap.scvi",
  group.by = c("conditions", "scvi_clusters"),
  combine = FALSE, label.size = 2, shuffle = T
) 

wrap_plots(c(p1, p2, p3, p4, p5), ncol = 5, byrow = F) & NoLegend()

p1 <- DimPlot(
  seurat_object,
  reduction = "umap.cca",
  group.by = "conditions",
  shuffle = T
) & NoLegend()
p2 <- DimPlot(
  seurat_object,
  reduction = "umap.rpca",
  group.by = "conditions",
  shuffle = T
) & NoLegend()
p3 <- DimPlot(
  seurat_object,
  reduction = "umap.harmony",
  group.by = "conditions",
  shuffle = T
) & NoLegend()
p4 <- DimPlot(
  seurat_object,
  reduction = "umap.mnn",
  group.by = "conditions",
  shuffle = T
) & NoLegend()
p5 <- DimPlot(
  seurat_object,
  reduction = "umap.scvi",
  group.by = "conditions",
  shuffle = T
) & NoLegend()

library(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, nrow = 1) 
arranged_plots <- grid.arrange(p1, p2, p3, p4, p5, nrow = 1)
arranged_plots

#save corrected object--------------------------------------------------------
saveRDS(seurat_object,
        "~/directory/seurat_object.rds")

#Now convert to h5ad for benchmarker-------------------------------------------
library(sceasy)
library(reticulate)

#join layers before assay conversion 
seurat_object <- readRDS("~/directory/seurat_object.rds")
seurat_object <- JoinLayers(seurat_object)

#convert V5 to V3 assay
seurat_object[["RNA"]] <- as(object = seurat_object[["RNA5"]], Class = "Assay")
DefaultAssay(seurat_object) <- "RNA"
unique(seurat_object@reductions)

#only preserve reductions that you would like to test
seurat_object <- DietSeurat(seurat_object, 
                                dimreducs = c("integrated.cca", "integrated.rpca", "harmony",
                                              "integrated.mnn", "integrated.scvi",
                                              "umap.cca", "umap.rpca", "umap.harmony",
                                              "umap.mnn", "umap.scvi"),
                                layers = c("counts", "data"))

sceasy::convertFormat(seurat_object, from="seurat", to="anndata",
                      outFile= '/Users/hauke/directory/seurat_object.h5ad')