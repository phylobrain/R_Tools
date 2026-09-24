# Automated annotations with Garnett

# Load packages
library(org.Gg.eg.db) # db for chicken
library(garnett)
library(ggplot2)

# LOAD THE DATA

# Seurat object to CDS
cds <- as.CellDataSet(seurat_object, 
                      assay = "SCT") # or RNA if not SCTransformed, expressionFamily = "negbinomial.size()" by default, recommended for most type of data

# Generate size factors for normalization later
cds <- estimateSizeFactors(cds)

# BUILD AND TRAIN THE CLASSIFIER 

# Construct the marker file and check the markers
marker_file_path <- file.choose()

marker_check <- check_markers(cds, marker_file_path,
                              db=org.Gg.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")

plot_markers(marker_check)

# Build and train the classifier
set.seed(260)
classifier <- train_cell_classifier(cds = cds,
                                    marker_file = marker_file_path,
                                    db=org.Gg.eg.db,
                                    cds_gene_id_type = "SYMBOL",
                                    num_unknown = 500,                          # 500 by default, outgroup cells for classification
                                    marker_file_gene_id_type = "SYMBOL")        # multinomial elastic net regression classification

# View the classification genes
feature_genes <- get_feature_genes(classifier,
                                   node = "root",                               # root for top node or cell type names for others
                                   db = org.Gg.eg.db)
head(feature_genes)

# View references 
get_classifier_references(classifier)

# CLASSIFY THE CELLS

# Classify the cells
cds <- classify_cells(cds, classifier,
                      db = org.Gg.eg.db,
                      cluster_extend = TRUE,
                      cds_gene_id_type = "SYMBOL")

# Look for the classifications
head(cds@phenoData@data$cell_type)
table(cds@phenoData@data$cell_type)

# Transfer to Seurat
seurat_object@meta.data$garnett_cell_type <- cds@phenoData@data$cell_type
DimPlot(seurat_object, group.by="garnett_cell_type")
