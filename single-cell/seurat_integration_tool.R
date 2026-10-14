# Seurat anchor based integration or classical integration
# Source: https://satijalab.org/seurat/articles/integration_rpca#performing-integration-on-datasets-normalized-with-sctransform

# Load packages
library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(dplyr)


# Import datasets
seurat_object1 <- readRDS("seurat_object1.rds")
seurat_object2 <- readRDS("seurat_object2.rds")

# Change default assay
DefaultAssay(seurat_object1) <- "RNA"
DefaultAssay(seurat_object2) <- "RNA"

# List the datasets that are going to be integrated
seurat_object_list <- list(seurat_object1, seurat_object2)


# SCTransform the objects of the list
seurat_object_list <- lapply(X = seurat_object_list, FUN = SCTransform, method = "glmGamPoi")

# Selection of the integration features
features <- SelectIntegrationFeatures(object.list = seurat_object_list)

# Run PCA on each dataset with the newly selected features
seurat_object_list <- PrepSCTIntegration(object.list = seurat_object_list, anchor.features = features)
seurat_object_list <- lapply(X = seurat_object_list, FUN = RunPCA, features = features)

# Identify the anchors
anchors <- FindIntegrationAnchors(object.list = seurat_object_list, 
                                  normalization.method = "SCT",
                                  anchor.features = features, 
                                  dims = 1:50,
                                  reduction = "cca") # Default; others include rpca...

# Integrate data
integration_object <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:50)

# Normal preprocessing, without repeating SCTransformation
integration_object <- RunPCA(integration_object, verbose = FALSE)
ElbowPlot(integration_object, ndims=50)
DimHeatmap(integration_object, dims = 1:20, cells = 500, balanced = TRUE, fast = F)

integration_object <- FindNeighbors(integration_object, reduction = "pca", dims = 1:15)
integration_object <- FindClusters(integration_object, resolution = c(0.5,1,1.5,2))
integration_object <- RunUMAP(integration_object, reduction = "pca", dims = 1:15)
