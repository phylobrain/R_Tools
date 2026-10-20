# LABEL TRANSFER TOOL
# Original source: https://satijalab.org/seurat/articles/integration_mapping

reference <- readRDS("object1.rds") # Annotated object to map the query dataset
query <- readRDS("object2.rds")

Idents(reference) <- "Cell_Identity" # variable with annotations for reference
Idents(query) <- "Cell_Identity" # variable with annotations for query

# data has been previously normalized or SCTransformed

# Find transfer anchors
anchors <- FindTransferAnchors(reference = reference, query = query, 
                               dims = 1:50,
                               reference.reduction = "pca",
                               normalization.method = "SCT")

# Transfer data and predict the new labels
predictions <- TransferData(anchorset = anchors, 
                            refdata = reference$Cell_Identity, 
                            dims = 1:50)

query <- AddMetaData(query, metadata = predictions) # add the predicted data to the query

DimPlot(query, group.by="predicted.id", pt.size = 1, label = T) # visualize predicted labels
table(query@meta.data$Cell_Identity, query@meta.data$predicted.id)

VlnPlot(query, c("LHX9", "EBF1"), group.by = "predicted.id") # view gene expression by predicted cell type

# Map query to reference UMAP structure
reference <- RunUMAP(reference, dims = 1:15, reduction = "pca", return.model = TRUE) # save the model
query <- MapQuery(anchorset = anchors, reference = reference, query = query,
                 refdata = list(celltype = "annotation"), # celltype = annotation variable name
                 reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(reference, reduction = "umap", group.by = "annotation", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE)+ NoLegend(ggtitle("Query transferred labels"))
p1 + p2


# Evaluation of predictions
DefaultAssay(query) <- "prediction.score.celltype"

FeaturePlot(query, features = unique(query$predicted.celltype),
            reduction = "ref.umap", cols = c("lightgrey", "darkred"), ncol=3)


                                                    