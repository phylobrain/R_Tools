# Annotations and visualization ready UMAP

# Change Idents
Idents(gg.integrated.neurons.E10) <- "SCT_snn_res.1.5"
gg.integrated.neurons.idents10 <- RenameIdents(gg.integrated.neurons.E10,
                                               c("0" = "",           
                                                 "1" = "",
                                                 "2" = "",        
                                                 "3" = "",       
                                                 "4" = "",
                                                 "5" = "",
                                                 "6" = "",     
                                                 "7" = "",
                                                 "8" = " ",
                                                 "9" = "",
                                                 "10" = "",
                                                 "11" = "",     
                                                 "12" = "",
                                                 "13" = "",     
                                                 "14" = "",               
                                                 "15" = "",
                                                 "16" = "",
                                                 "17" = "", 
                                                 "18" = "",            
                                                 "19" = "",
                                                 "20" = "",
                                                 "21" = "",      
                                                 "22" = "",
                                                 "23" = "",
                                                 "24" = "",
                                                 "25" = "",     
                                                 "26" = "", 
                                                 "27" = "",
                                                 "28" = "",         
                                                 "29" = "",               
                                                 "30" = "",    
                                                 "31" = "",
                                                 "32" = "",            
                                                 "33" = "",
                                                 "34" = "",         
                                                 "35" = "",               
                                                 "36" = ""))            

# Visualization ready UMAP by colors
umap_df <- as.data.frame(Embeddings(gg.integrated.neurons.idents10, "umap"))
umap_df$cluster <- Idents(gg.integrated.neurons.idents10) # get the coordinates

centers <- umap_df %>%      # set the center for each cluster
  group_by(cluster) %>%
  summarise(
    UMAP_1 = mean(umap_1),
    UMAP_2 = mean(umap_2)
  )

cluster_colors <- c(
  "" = "#AB3F6A",        
  "" = "#6A6A6A",
  "" = "#B494D5",        
  "" = "#B0A9D3",      
  "" = "#610070",
  "" = "#90301B",
  "" = "#A3C395",      
  "" = "#376491",
  "" = "#7CADD2",
  "" = "#3D1778",
  "" = "#454545",
  "" = "#81B357",     
  "" = "#6851A5",
  "" = "#CBE59E",     
  "" = "#AE123A",               
  "" = "#F36D60",
  "" = "#ED5F54",
  "" = "#46924F", 
  "" = "#24693D",            
  "" = "#398347",
  "" = "#7665AE",
  "" = "#AB88BA",      
  "" = "#C42F40",
  "" = "#858585",
  "" = "#D54046",
  "" = "#6FAA66",    
  "" = "#AADAD4", 
  "" = "#93AD5C",
  "" = "#B84585",          
  "" = "#E9982E",               
  "" = "#5E9B7C",   
  "" = "#6E0278",
  "" = "#6E9F6D",            
  "" = "#358AAA",
  "" = "#E8A3B8",           
  "" = "#54BDC2",              
  "" = "#C8A7E1",    
  "" = "#E3C181",
  "" = "#533492",
  "" = "#A8572B",
  "" = "#F5C26E",           
  "" = "#8576B7",
  "" = "#7664AE",
  "" = "#86CAE1"
)

DimPlot(
  gg.integrated.neurons.idents10,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.5,
  cols = cluster_colors
) +
  geom_label_repel(
    data = centers,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = cluster,
      fill = cluster,
    ),
    color = "white",
    size = 2.5,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "black",
    segment.size = 0.4
  ) +
  scale_fill_manual(values = cluster_colors) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

DimPlot(
  gg.integrated.neurons.idents10,
  reduction = "umap",
  label = FALSE,
  pt.size = 0.5,
  cols = cluster_colors
) +
  geom_label_repel(
    data = centers,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = cluster,
      color = cluster
    ),
    fill = scales::alpha("white",0.8),
    size = 2.5,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = cluster_colors) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

# Add new idents as a new column
gg.integrated.neurons.E10$Cell_Identity <- Idents(gg.integrated.neurons.idents10)
DimPlot(gg.integrated.neurons.E10, group.by = "Cell_Identity", label = T)
