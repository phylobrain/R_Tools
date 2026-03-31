# Load packages
library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(dplyr)

# Import datasets
seurat_object1 <- readRDS("seurat_object1.rds")
seurat_object2 <- readRDS("seurat_object2.rds")


# Seurat anchor based integration

# List of the Seurat objects
seurat_object_list <- list(seurat_object1, seurat_object2)

# Selection of the integration features
features <- SelectIntegrationFeatures(object.list = seurat_object_list)

# Find integration anchors
seurat_object_list <- PrepSCTIntegration(object.list = seurat_object_list,
                                         anchor.features = features) # Only when needed because SCT transformation in separate objects

anchors <- FindIntegrationAnchors(object.list = seurat_object_list, 
                                  anchor.features = features,
                                  scale = TRUE,
                                  normalization.method = "SCT",
                                  reduction = "cca") # Default, other options include "rpca", etc.

# Integration
integrated_object <- IntegrateData(anchorset = anchors)


# Necessary batch effect analysis via agrupation given the origin of samples, species... or other confounding variables

# For dataset preprocessing and clusterization change the default assay to SCT or RNA
DefaultAssay(integrated_object) <- "SCT" # or RNA

# The new dataset must be preprocessed, clusterd and annotated as the individual objects

