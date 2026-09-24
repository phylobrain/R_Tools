# Install monocle3
remotes::install_github("bnprks/BPCells/r")
devtools::install_github('cole-trapnell-lab/monocle3')

# Load packages
library(monocle3)
library(Seurat)
library(SeuratWrappers)


# Load data
seurat_object <- readRDS("")

# Seurat to CDS
cds <- as.cell_data_set(seurat_object, 
                      assay = "SCT") # or RNA if not SCTransformed, expressionFamily = "negbinomial.size()" by default, recommended for most type of data


# Plot
plot_cells(cds, color_cells_by = "") # color by any column in colData(cds), same UMAP as in Seurat
plot_cells(cds, genes = "") # as FeaturePlot()


# OPTIONAL - Remove batch effect (coming from Seurat probably it is corrected)
cds <- align_cds(cds, num_dim = 50, alignment_group = "") # Default preprocess_method = "PCA"
cds <- reduce_dimension(cds)


# Cluster cells
cds <- cluster_cells(cds) # by default uses Levine community detection algorithm, essential to get cluster partitions

# Finding marker genes by each cluster
marker_test_res <- top_markers(cds,
                               group_cells_by = "") # Based on the differential testing

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2) # Select top n markers

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    max.size=3)

# Annotations
colData(cds)$assigned_cell_type <- as.character(cluster(cds))

colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$assigned_cell_type,
                                                 "1"="",
                                                 "2"="",
                                                 "3"="",
                                                 "4"="",
                                                 "5"="",
                                                 "6"="",
                                                 "7"="",
                                                 "8"="")

plot_cells(cds, group_cells_by="cluster", color_cells_by="assigned_cell_type")

# Isolate cells of interest for further analysis
cds_subset <- choose_cells(cds)

# DE
pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="knn", cores=8) # identify DE genes in the sybclusters
pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=1e-3) # group DE that have a similar pattern of expression
plot_cells(cds_subset, genes=gene_module_df, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE) # plot the modules


# CONSTRUCTING TRAJECTORIES

# Prepare the data
cds <- as.cell_data_set(seurat_object)

cds <- cluster_cells(cds)

# Learn the graph
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

plot_cells(cds,
           color_cells_by = "", # time to select where the trajectory starts
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

# Order the cells
cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5) # color by pseudotime

# Helper function to identify the root principal points
get_earliest_principal_node <- function(cds, time_bin=""){ # earliest time of your data
  cell_ids <- which(colData(cds)[, ""] == time_bin) # variable that measures developmental time in your data
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

# Subset by branch
cds_sub <- choose_graph_segments(cds)

# 3D trajectory graphs
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")

# Saving and loading monocle objects
save_monocle_objects(cds = cds, directory_path = "cds_object")
cds <- load_monocle_objects(directory_path = "cds_object")





