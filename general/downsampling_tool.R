# Load packages
library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(dplyr)

# Import data
small_obj <- readRDS("")
big_obj <- readRDS("")

# Stratified downsampling without replacement of the big object to fit the small object size

# Set the proportions to mantain
prop <- prop.table(table(small_obj@meta.data$cell_type)) # cell_type or group to mantain props
n_target <- length(colnames(small_obj))
n_by_group <- round(n_target * prop)

set.seed(123) # random seed for reproducibility
downsampling <- unlist(lapply(names(n_by_group), function(x){                   # selection of cells without replacement
  sample(
    WhichCells(big_obj, expression = annotation == x),
    size = n_by_group[x]
  )
}))

downsampled_obj <- subset(big_obj, cells = downsampling)