if (!exists("snakemake")) {
  stop("Run this script through Snakemake.")
}

results_dir <- snakemake@output[[1]]
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
setwd(results_dir)

data <- readRDS(snakemake@input[["annotation_data"]])
data_con <- subset(data, group %in% grep("_con$", group, value = TRUE))
neurons <- readRDS(snakemake@input[["neurons"]])
neurons_con <- subset(neurons, group %in% grep("_con$", group, value = TRUE))

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
# Cell-class composition and t-SNE
library(Seurat)
library(dplyr)
library(ggplot2)




#data <- readRDS("annotation_data.rds")
#data_con <-  data %>%
  #subset(group %in% grep("_con$", group, value = TRUE))

name_dict <- c(
  "Acquisition_con" = "td+ey-",
  "Retrieval_con" = "td-ey+",
  "Overlapping_con" = "td+ey+",
  "Other_con" = "td-ey-"
)
data_con@meta.data$group <- recode(data_con@meta.data$group, !!!name_dict)



meta_data <- data_con@meta.data

cell_counts <- meta_data %>%
  group_by(group, class) %>%
  summarise(count = n(), .groups = 'drop')


total_counts <- cell_counts %>%
  group_by(group) %>%
  summarise(total = sum(count), .groups = 'drop')


proportion_data <- cell_counts %>%
  left_join(total_counts, by = "group") %>%
  mutate(proportion = count / total)


#proportion
custom_colors <- c("#A5D1B0", "#CE8A8D", "#FFF7C1", "#E0F3FF", 
                   "#ADD3F4", "#F7C9CF", "#FEE4E8", "#7CA3B8",
                   "#BFB8D6", "#FCCB8E", "#FFA07A")
proportion_data$subclass <- factor(proportion_data$class, 
                                   levels = c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", 
                                              "Endothelial cells", "Fibroblast", "Monocyte", "Macrophage",
                                              "T cells", "B cell", "NK"))
proportion_data$group <- factor(proportion_data$group, 
                                levels = c("td-ey-", 
                                           "td+ey-", 
                                           "td-ey+", 
                                           "td+ey+"))
plotproportion <- ggplot(proportion_data, aes(x = group, y = proportion, fill = class)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(x = "Group", y = "Proportion", fill = "class") +
  theme(axis.text.x = element_text(size = 40, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40)
  )

ggsave(
  filename = paste0("E:/master_graduation/tsne_proportion/proportion.png"),
  plot = plotproportion,      
  device = "png",              
  width = 3840,           
  height = 2560,         
  units = "px",                 
  dpi = 300,
  limitsize = FALSE
)


#tsne
tsne_data <- as.data.frame(data_con@reductions$tsne@cell.embeddings)
tsne_data$cluster <- data_con@meta.data$class 
tsne_data$cluster <- factor(tsne_data$cluster, 
                            levels = c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", 
                                       "Endothelial cells", "Fibroblast", "Monocyte", "Macrophage",
                                       "T cells", "B cell", "NK"))
ptsne <- ggplot(tsne_data, aes(x = tSNE_1, y = tSNE_2, color = cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  labs(x = "t-SNE 1", y = "t-SNE 2", color = "class") +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 40),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40)
  )

ggsave(
  filename = paste0("E:/master_graduation/tsne_proportion/tsne.png"),
  plot = ptsne,      
  device = "png",              
  width = 3840,           
  height = 2560,         
  units = "px",                 
  dpi = 300,
  limitsize = FALSE
)


# Neuronal-subclass composition and t-SNE
library(Seurat)
library(dplyr)
library(ggplot2)



meta_data <- neurons_con@meta.data

cell_counts <- meta_data %>%
  group_by(group, subclass) %>%
  summarise(count = n(), .groups = 'drop')


total_counts <- cell_counts %>%
  group_by(group) %>%
  summarise(total = sum(count), .groups = 'drop')


proportion_data <- cell_counts %>%
  left_join(total_counts, by = "group") %>%
  mutate(proportion = count / total)

group_mapping <- data.frame(  
  actual = c("Other_con", "Acquisition_con", "Retrieval_con", "Overlapping_con"),  
  display = c("td-ey-", "td+ey-", "td-ey+", "td+ey+")  
)  

proportion_data <- proportion_data %>%  
  left_join(group_mapping, by = c("group" = "actual")) %>%  
  mutate(group = display)

#proportion
custom_colors <- c("#A5D1B0", "#CE8A8D", "#FFF7C1", "#E0F3FF", 
                   "#ADD3F4", "#F7C9CF", "#FEE4E8", "#7CA3B8", "#BFB8D6")
proportion_data$subclass <- factor(proportion_data$subclass, 
                                   levels = c("C0", "C1", "C2", "C3", "C4", 
                                              "C5", "C6", "C7"))
proportion_data$group <- factor(proportion_data$group, 
                                levels = c("td-ey-", 
                                           "td+ey-", 
                                           "td-ey+", 
                                           "td+ey+"))
plotproportion <- ggplot(proportion_data, aes(x = group, y = proportion, fill = subclass)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(x = "Group", y = "Proportion", fill = "Subclass") +
  theme(axis.text.x = element_text(size = 40, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40)
  )

ggsave(
  filename = paste0("E:/master_graduation/tsne_proportion/proportion_subclass.png"),
  plot = plotproportion,      
  device = "png",              
  width = 3840,           
  height = 2560,         
  units = "px",                 
  dpi = 300,
  limitsize = FALSE
)


#tsne
tsne_data <- as.data.frame(neurons_con@reductions$tsne@cell.embeddings)
tsne_data$cluster <- neurons_con@meta.data$subclass 
tsne_data$cluster <- factor(tsne_data$cluster, 
                            levels = c("C0", "C1", "C2", "C3", "C4", 
                                       "C5", "C6", "C7"))
ptsne <- ggplot(tsne_data, aes(x = tSNE_1, y = tSNE_2, color = cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  labs(x = "t-SNE 1", y = "t-SNE 2", color = "Subclass") +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 40),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40)
  )

ggsave(
  filename = paste0("E:/master_graduation/tsne_proportion/tsne_subclass.png"),
  plot = ptsne,      
  device = "png",              
  width = 3840,           
  height = 2560,         
  units = "px",                 
  dpi = 300,
  limitsize = FALSE
)
