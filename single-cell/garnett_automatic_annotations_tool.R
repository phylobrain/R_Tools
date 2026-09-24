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
marker_file_path <- system.file("extdata", "markers.txt",
                                package = "garnett")

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

# Plot the classification results
qplot(umap_1, umap_2, color = cell_type, data = pData(cds))    
qplot(umap_1, umap_2, color = cluster_ext_type, data = pData(cds))
