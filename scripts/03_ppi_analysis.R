if (!exists("snakemake")) {
  stop("Run this script through Snakemake.")
}

results_dir <- snakemake@output[[1]]
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
setwd(results_dir)

map_output_path <- function(path) {
  if (grepl("^E:/master_graduation", path)) {
    return(file.path(results_dir, sub("^E:/master_graduation/?", "", path)))
  }
  path
}

ggsave <- function(filename, ...) {
  filename <- map_output_path(filename)
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = filename, ...)
}
library(Seurat)
library(dplyr)
library(tibble)

input_deg_dir <- file.path(snakemake@input[["neuronal_results"]], "differential_expression")
read_deg_table <- function(filename) {
  read.delim(file.path(input_deg_dir, filename), row.names = 1, check.names = FALSE)
}

acquisition_filtered <- read_deg_table("acquisition_vs_other.tsv")
retrieval_filtered <- read_deg_table("retrieval_vs_other.tsv")
overlapping_filtered <- read_deg_table("overlapping_vs_other.tsv")
dir.create(file.path(results_dir, "PPI"), recursive = TRUE, showWarnings = FALSE)

# PPI networks
library(STRINGdb)  
library(igraph)  
library(ggraph)  
library(ggplot2)  


string_db <- STRINGdb$new(version = "12.0", species = 10090, score_threshold = 400, input_directory = "")  


acquisition_genes <- data.frame(gene = rownames(acquisition_filtered),  
                                avg_log2FC = acquisition_filtered$avg_log2FC)  

retrieval_genes <- data.frame(gene = rownames(retrieval_filtered),  
                              avg_log2FC = retrieval_filtered$avg_log2FC)  

overlapping_genes <- data.frame(gene = rownames(overlapping_filtered),  
                                avg_log2FC = overlapping_filtered$avg_log2FC)  


acq_mapped <- string_db$map(acquisition_genes, "gene", removeUnmappedRows = TRUE)  
ret_mapped <- string_db$map(retrieval_genes, "gene", removeUnmappedRows = TRUE)  
ovl_mapped <- string_db$map(overlapping_genes, "gene", removeUnmappedRows = TRUE)  


get_largest_connected_component <- function(graph) {  

  comps <- components(graph)  

  largest_comp <- which.max(comps$csize)  
  return(induced_subgraph(graph, which(comps$membership == largest_comp)))  
}  


get_ppi_graph <- function(mapped_data, output_path) {  

  network <- string_db$get_interactions(mapped_data$STRING_id)  
  

  igraph_network <- graph_from_data_frame(network, directed = FALSE)  
  

  igraph_network <- get_largest_connected_component(igraph_network)  
  

  V(igraph_network)$avg_log2FC <- mapped_data$avg_log2FC[match(V(igraph_network)$name, mapped_data$STRING_id)]  
  V(igraph_network)$label <- mapped_data$gene[match(V(igraph_network)$name, mapped_data$STRING_id)]  
  

  V(igraph_network)$degree <- degree(igraph_network)  
  

  V(igraph_network)$group <- ifelse(V(igraph_network)$avg_log2FC > 0, "up", "down")  
  
  hub_nodes <- order(V(igraph_network)$degree, decreasing = TRUE)[1:10]
  hub_genes <- V(igraph_network)$label[hub_nodes]
  V(igraph_network)$label_to_show <- ifelse(V(igraph_network)$label %in% hub_genes, V(igraph_network)$label, "")
  

  groups_present <- unique(V(igraph_network)$group)
  

  node_colors <- c("up" = "#C30078", "down" = "#4370B4")
  scale_labels <- c()
  

  if ("up" %in% groups_present & "down" %in% groups_present) {  
    scale_labels <- c("Down", "Up")  
  } else if ("up" %in% groups_present) {  
    node_colors <- node_colors["up"]  
    scale_labels <- c("Up")  
  } else if ("down" %in% groups_present) {  
    node_colors <- node_colors["down"]  
    scale_labels <- c("Down")  
  }  
  

  plot <- ggraph(igraph_network, layout = "kk") +  

    geom_edge_link(edge_colour = "gray", edge_width = 0.5, alpha = 0.6) +   

    geom_node_point(aes(fill = group, size = degree), shape = 21) +   

    #geom_node_text(aes(label = label_to_show), repel = TRUE, size = 11.5) +  

    scale_fill_manual(  
      values = node_colors,
      name = "Expression",
      labels = scale_labels
    ) +   

    scale_size_continuous(range = c(5, 20), name = "Degree") +  
    theme_void() +
    theme(  
      legend.text = element_text(size = 80),
      legend.title = element_text(size = 80, margin = margin(b = 20)),
      #legend.key.size = unit(10, "cm")
      legend.margin = margin(b = 20)
    ) 
  

  ggsave(  
    filename = output_path,  
    plot = plot,  
    device = "png",              
    width = 3840 * 2,
    height = 2560 * 2,           
    units = "px",
    dpi = 300
  )  
}  


get_ppi_graph(  
  mapped_data = acq_mapped,  
  output_path = "E:/master_graduation/PPI/acq.png"  
)  


get_ppi_graph(  
  mapped_data = ret_mapped,  
  output_path = "E:/master_graduation/PPI/ret.png"  
)  


get_ppi_graph(  
  mapped_data = ovl_mapped,  
  output_path = "E:/master_graduation/PPI/ovl.png"  
)

# PPI betweenness visualizations
library(STRINGdb)  
library(igraph)  
library(ggraph)  
library(ggplot2)  


string_db <- STRINGdb$new(version = "12.0", species = 10090, score_threshold = 400, input_directory = "")  


get_largest_connected_component <- function(graph) {  

  comps <- components(graph)  

  largest_comp <- which.max(comps$csize)  
  return(induced_subgraph(graph, which(comps$membership == largest_comp)))  
}  


get_ppi_graph_with_gradients <- function(mapped_data, output_path) {  

  network <- string_db$get_interactions(mapped_data$STRING_id)  
  

  igraph_network <- graph_from_data_frame(network, directed = FALSE)  
  

  igraph_network <- get_largest_connected_component(igraph_network)  
  

  V(igraph_network)$avg_log2FC <- mapped_data$avg_log2FC[match(V(igraph_network)$name, mapped_data$STRING_id)]  
  V(igraph_network)$label <- mapped_data$gene[match(V(igraph_network)$name, mapped_data$STRING_id)]  
  

  V(igraph_network)$betweenness <- betweenness(igraph_network, normalized = TRUE)
  V(igraph_network)$degree <- degree(igraph_network)
  

  hub_nodes <- order(V(igraph_network)$betweenness, decreasing = TRUE)[1:10]
  hub_genes <- V(igraph_network)$label[hub_nodes]
  V(igraph_network)$label_to_show <- ifelse(V(igraph_network)$label %in% hub_genes, V(igraph_network)$label, "")
  

  gradient_colors <- c("black", "purple", "orange", "wheat")  
  

  plot <- ggraph(igraph_network, layout = "kk") +  

    geom_edge_link(edge_colour = "gray", edge_width = 0.5, alpha = 0.6) +   

    geom_node_point(aes(fill = betweenness, size = degree), shape = 21) +   

    #geom_node_text(aes(label = label_to_show), color = "black", repel = TRUE, size = 11.5) +  

    scale_fill_gradientn(  
      colors = gradient_colors,  
      name = "Betweenness"
    ) +  

    scale_size_continuous(range = c(5, 20), name = "Degree") +  

    theme_void() +  
    theme(  
      legend.position = "right",
      legend.title = element_text(size = 80, margin = margin(b = 20)),  
      legend.text = element_text(size = 80),
      legend.key.height = unit(2, "cm"),
      legend.margin = margin(b = 20)
    )  
  

  ggsave(  
    filename = output_path,  
    plot = plot,  
    device = "png",              
    width = 3840 * 2,
    height = 2560 * 2,           
    units = "px",
    dpi = 300
  )  
}  


get_ppi_graph_with_gradients(  
  mapped_data = acq_mapped,  
  output_path = "E:/master_graduation/PPI/acq_hub_gradient.png"  
)  


get_ppi_graph_with_gradients(  
  mapped_data = ret_mapped,  
  output_path = "E:/master_graduation/PPI/ret_hub_gradient.png"  
)  


get_ppi_graph_with_gradients(  
  mapped_data = ovl_mapped,  
  output_path = "E:/master_graduation/PPI/ovl_hub_gradient.png"  
)

# PPI hub-gene tables
library(STRINGdb)  
library(igraph)  
library(openxlsx)


string_db <- STRINGdb$new(version = "12.0", species = 10090, score_threshold = 400, input_directory = "")  


get_largest_connected_component <- function(graph) {  

  comps <- components(graph)  

  largest_comp <- which.max(comps$csize)  
  return(induced_subgraph(graph, which(comps$membership == largest_comp)))  
}  


extract_and_save_hub_genes <- function(mapped_data, xlsx_path) {  

  network <- string_db$get_interactions(mapped_data$STRING_id)  
  

  igraph_network <- graph_from_data_frame(network, directed = FALSE)  
  

  igraph_network <- get_largest_connected_component(igraph_network)  
  

  V(igraph_network)$avg_log2FC <- mapped_data$avg_log2FC[match(V(igraph_network)$name, mapped_data$STRING_id)]  
  V(igraph_network)$label <- mapped_data$gene[match(V(igraph_network)$name, mapped_data$STRING_id)]  
  

  V(igraph_network)$betweenness <- betweenness(igraph_network, normalized = TRUE)
  V(igraph_network)$degree <- degree(igraph_network)
  

  degree_top20 <- order(V(igraph_network)$degree, decreasing = TRUE)[1:30]  
  betweenness_top20 <- order(V(igraph_network)$betweenness, decreasing = TRUE)[1:30]  
  

  intersect_genes <- intersect(V(igraph_network)$label[degree_top20], V(igraph_network)$label[betweenness_top20])  
  hub_genes <- head(intersect_genes, 10)
  

  print("Hub genes (intersection of top 20 Degree and Betweenness):")  
  print(hub_genes)  
  

  hub_genes_df <- data.frame(Gene = hub_genes)
  write.xlsx(hub_genes_df, xlsx_path)
}  




extract_and_save_hub_genes(  
  mapped_data = acq_mapped,  
  xlsx_path = map_output_path("E:/master_graduation/PPI/acq_hub_genes.xlsx")
)  


extract_and_save_hub_genes(  
  mapped_data = ret_mapped,  
  xlsx_path = map_output_path("E:/master_graduation/PPI/ret_hub_genes.xlsx")  
)  


extract_and_save_hub_genes(  
  mapped_data = ovl_mapped,  
  xlsx_path = map_output_path("E:/master_graduation/PPI/ovl_hub_genes.xlsx")  
)
