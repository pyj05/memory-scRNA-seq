library(Seurat)
library(NeuronChat)
library(dplyr)
library(CellChat)
library(cowplot)
library(grid)
library(ComplexHeatmap)
library(RColorBrewer)

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
neurons <- readRDS(snakemake@input[["neurons"]])
neurons_con <- subset(neurons, group %in% grep("_con$", group, value = TRUE))

# NeuronChat analysis across experimental groups
subset_Acquisition_con <- subset(neurons_con, subset = group == 'Acquisition_con')
subset_Retrieval_con <- subset(neurons_con, subset = group == 'Retrieval_con')
subset_Overlapping_con <- subset(neurons_con, subset = group == 'Overlapping_con')
subset_Other_con <- subset(neurons_con, subset = group == 'Other_con')

set.seed(1234)
seurat_list <- list(
  'td-ey-' = subset_Other_con,
  'td+ey-' = subset_Acquisition_con,
  'td-ey+' = subset_Retrieval_con,
  'td+ey+' = subset_Overlapping_con
)
neuronchat_list <- lapply(names(seurat_list), function(group_name) {
  seurat_obj <- seurat_list[[group_name]]
  target_df <- GetAssayData(seurat_obj, layer = "data") 
  meta <- seurat_obj@meta.data
  colnames(target_df) <- rownames(meta)
  meta <- meta[!is.na(meta$subclass), ]
  
  neuronchat_obj <- createNeuronChat(
    as.matrix(target_df),  
    DB = "mouse",                   
    group.by = meta$subclass,       
    meta = meta                     
  )
  
  neuronchat_obj <- run_NeuronChat(neuronchat_obj, M = 100)
  return(neuronchat_obj)
})
names(neuronchat_list) <- names(seurat_list)


file_name <- map_output_path("E:/master_graduation/neuronchatplot/netVisual_circle_neuron.pdf")
dir.create(dirname(file_name), recursive = TRUE, showWarnings = FALSE)
pdf(file = file_name, width = 15.36, height = 8.64)  
par(mfrow = c(1, 4))
for (j in c(1, 2, 3, 4)) {
  net_aggregated_x <- net_aggregation(neuronchat_list[[j]]@net, method = 'weight')
  netVisual_circle_neuron(net_aggregated_x, title.name = names(neuronchat_list)[j], arrow.size = 0.5, margin = 0.4, edge.width.max = 10, vertex.label.cex = 3)
}
dev.off()



neuronchat_merge_list <- mergeNeuronChat(neuronchat_list, add.names = names(neuronchat_list))

my_colors <- c("#D5D9E5", "#F7A6AC", "#B2DBB9", "#EEF0A7")

p1 <- compareInteractions_Neuron(neuronchat_merge_list,measure = c("weight"),comparison = c(1,2,3,4),group=c(1,2,3,4),show.legend = F, color.use = my_colors, size.text = 60, width = 0.1, x.lab.rot = TRUE)
p2 <- compareInteractions_Neuron(neuronchat_merge_list,measure = c("count"),comparison = c(1,2,3,4),group=c(1,2,3,4),show.legend = F, color.use = my_colors, size.text = 60, width = 0.1, x.lab.rot = TRUE)
p1 <- p1 + scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(plot.margin = unit(c(0,5,0,0),units= "cm"))
p2 <- p2 + scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) 
ggsave(
  filename = "E:/master_graduation/neuronchatplot/compareInteractions_Neuron.pdf",
  plot = p1+p2,      
  device = "pdf",              
  width = 3840 * 2,           
  height = 2560 * 2,         
  units = "px",                 
  dpi = 300,
  limitsize = FALSE
)


g1 <- rankNet_Neuron(neuronchat_merge_list,mode='comparison',measure = c("count"),comparison = c(1,2,3,4),do.stat = F,tol = 0.1,stacked = F,font.size = 30, color.use = my_colors)
g2 <- rankNet_Neuron(neuronchat_merge_list,mode='comparison',measure = c("weight"),comparison = c(1,2,3,4),do.stat = F,tol = 0.1,stacked = F,font.size = 30, color.use = my_colors) 
g1 <- g1 + theme(legend.text = element_text(size = 40),
                 plot.margin = unit(c(0,5,0,0),units= "cm"),
                 legend.key.size = unit(2, "lines"),
                 axis.title.x = element_text(size = 40))
g2 <- g2 + theme(legend.text = element_text(size = 40),
                 legend.key.size = unit(2, "lines"),
                 axis.title.x = element_text(size = 40))
ggsave(
  filename = "E:/master_graduation/neuronchatplot/rankNet_Neuron.pdf",
  plot = g1+g2,      
  device = "pdf",              
  width = 3840 * 3,           
  height = 2560 * 3,         
  units = "px",                 
  dpi = 300,
  limitsize = FALSE
)



all_cell_types <- paste0("C", 0:7)


expand_matrix <- function(mat, all_cell_types) {

  full_mat <- matrix(0, nrow = length(all_cell_types), ncol = length(all_cell_types))
  rownames(full_mat) <- all_cell_types
  colnames(full_mat) <- all_cell_types
  

  common_rows <- intersect(rownames(mat), all_cell_types)
  common_cols <- intersect(colnames(mat), all_cell_types)
  full_mat[common_rows, common_cols] <- mat[common_rows, common_cols]
  
  return(full_mat)
}


acquisition_con_list <- neuronchat_merge_list@net$`td+ey-`
for (i in 1:length(acquisition_con_list)) {
  acquisition_con_list[[i]] <- expand_matrix(acquisition_con_list[[i]], all_cell_types)
}
neuronchat_merge_list@net$`td+ey-` <- acquisition_con_list


retrieval_con_list <- neuronchat_merge_list@net$`td-ey+`
for (i in 1:length(retrieval_con_list)) {
  retrieval_con_list[[i]] <- expand_matrix(retrieval_con_list[[i]], all_cell_types)
}
neuronchat_merge_list@net$`td-ey+` <- retrieval_con_list


overlapping_con_list <- neuronchat_merge_list@net$`td+ey+`
for (i in 1:length(overlapping_con_list)) {
  overlapping_con_list[[i]] <- expand_matrix(overlapping_con_list[[i]], all_cell_types)
}
neuronchat_merge_list@net$`td+ey+` <- overlapping_con_list


other_con_list <- neuronchat_merge_list@net$`td-ey-`
for (i in 1:length(other_con_list)) {
  other_con_list[[i]] <- expand_matrix(other_con_list[[i]], all_cell_types)
}
neuronchat_merge_list@net$`td-ey-` <- other_con_list


#compute functional similarity 
neuronchat_merge_list <- computeNetSimilarityPairwise_Neuron(neuronchat_merge_list, slot.name = "net", type = "functional",comparison = c(1,2,3,4))
#Compute signaling network similarity for datasets 
#manifold learning
neuronchat_merge_list <- netEmbedding(neuronchat_merge_list,slot.name = "net_analysis", type = "functional",comparison = c(1,2,3,4))
#Manifold learning of the signaling networks for datasets 
#clustering on interactions 
neuronchat_merge_list <- netClustering(neuronchat_merge_list, slot.name = "net_analysis", type = "functional",comparison = c(1,2,3,4),k = 5)
#Classification learning of the signaling networks for datasets 

#add "edge.stroke" parameter
netVisual_embeddingPairwise_Neuron <- function(object, slot.name = "net_analysis", type = c("functional","structural"), comparison = NULL, color.use = NULL, point.shape = NULL, pathway.labeled = NULL, top.label = 1, pathway.remove = NULL, pathway.remove.show = TRUE, dot.size = c(2, 6), label.size = 2.5, dot.alpha = 0.5,
                                               xlabel = "Dim 1", ylabel = "Dim 2", title = NULL,do.label = T, show.legend = T, show.axes = T, edge.stroke = 1) {
  type <- match.arg(type)
  if (is.null(comparison)) {
    comparison <- 1:length(unique(object@meta$datasets))
  }
  cat("2D visualization of signaling networks from datasets", as.character(comparison), '\n')
  comparison.name <- paste(comparison, collapse = "-")
  
  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  clusters <- methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]]
  object.names <- names(methods::slot(object, 'net'))[comparison]
  prob <- list()
  for (i in 1:length(comparison)) {
    object.net <- methods::slot(object, 'net')[[comparison[i]]]
    prob[[i]] = simplify2array(object.net)
  }
  names(prob) <-object.names
  
  if (is.null(point.shape)) {
    point.shape <- 1:15
  }
  
  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
    pathway.remove <- rownames(similarity)[which(colSums(similarity) == 1)]
    pathway.remove <- sub("--.*", "", pathway.remove)
  }
  
  if (length(pathway.remove) > 0) {
    for (i in 1:length(prob)) {
      probi <- prob[[i]]
      pathway.remove.idx <- which(paste0(dimnames(probi)[[3]],"--",object.names[i]) %in% pathway.remove)
      if (length(pathway.remove.idx) > 0) {
        probi <- probi[ , , -pathway.remove.idx]
      }
      prob[[i]] <- probi
    }
  }
  prob_sum.each <- list()
  signalingAll <- c()
  for (i in 1:length(prob)) {
    probi <- prob[[i]]
    prob_sum.each[[i]] <- apply(probi, 3, sum)
    signalingAll <- c(signalingAll, paste0(names(prob_sum.each[[i]]),"--",object.names[i]))
  }
  prob_sum <- unlist(prob_sum.each)
  names(prob_sum) <- signalingAll
  prob_sum <- prob_sum[rownames(Y)]
  
  group <- sub(".*--", "", names(prob_sum))
  labels = sub("--.*", "", names(prob_sum))
  
  df <- data.frame(x = Y[,1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum),
                   labels = as.character(labels), clusters = as.factor(clusters), group = factor(group, levels = unique(group)))
  
  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(clusters)))
  }
  gg <- ggplot(data = df, aes(x, y)) +
    geom_point(aes(size = Commun.Prob., fill = clusters, colour = clusters, shape = group),
               stroke = edge.stroke) +
    CellChat_theme_opts() +
    theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in"))+
    guides(colour = guide_legend(override.aes = list(size = 3)))+
    labs(title = title, x = xlabel, y = ylabel) +
    scale_size_continuous(limits = c(0,1), range = dot.size, breaks = c(0.1,0.5,0.9)) +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
    theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25))
  gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, alpha = dot.alpha), drop = FALSE)
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE)
  gg <- gg + scale_shape_manual(values = point.shape[1:length(prob)])
  if (do.label) {
    gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, colour = clusters, alpha=group), size = label.size, show.legend = F,segment.size = 0.2, segment.alpha = 0.5) + scale_alpha_discrete(range = c(1, 0.6))
  }
  
  if (length(pathway.remove) > 0 & pathway.remove.show) {
    gg <- gg + annotate(geom = 'text', label =  paste("Isolate pathways: ", paste(pathway.remove, collapse = ', ')), x = -Inf, y = Inf, hjust = 0, vjust = 1, size = label.size,fontface="italic")
  }
  
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  
  if (!show.axes) {
    gg <- gg + theme_void()
  }
  gg
}


#embeddingPairwise
p3 <- netVisual_embeddingPairwise_Neuron(neuronchat_merge_list, slot.name = "net_analysis", type = "functional", label.size = 20,comparison=c(1,2,3,4),pathway.remove.show = FALSE,pathway.labeled = F, dot.size = c(10, 20), edge.stroke = 3)
p3 <- p3 + theme(legend.text = element_text(size = 40),
                 legend.title = element_text(size = 40),
                 legend.key.size = unit(2, "lines"),
                 axis.title.x = element_text(size = 40),
                 axis.title.y = element_text(size = 40))
ggsave(
  filename = "E:/master_graduation/neuronchatplot/embeddingPairwise.pdf",
  plot = p3,      
  device = "pdf",              
  width = 3840 * 3,           
  height = 2560 * 3,         
  units = "px",                 
  dpi = 300,
  limitsize = FALSE
)


#add "stroke.size"parameter
netVisual_embeddingPairwiseZoomIn_Neuron <- function(object, slot.name = "net_analysis", type = c("functional","structural"), comparison = NULL, color.use = NULL, nCol = 1, point.shape = NULL, pathway.remove = NULL, dot.size = c(2, 6), label.size = 2.8, dot.alpha = 0.5,
                                                     stroke.size = 0.5,
                                                     xlabel = NULL, ylabel = NULL, do.label = T, show.legend = F, show.axes = T) {
  type <- match.arg(type)
  if (is.null(comparison)) {
    comparison <- 1:length(unique(object@meta$datasets))
  }
  cat("2D visualization of signaling networks from datasets", as.character(comparison), '\n')
  comparison.name <- paste(comparison, collapse = "-")
  
  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  clusters <- methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]]
  object.names <- names(methods::slot(object, 'net'))[comparison]
  prob <- list()
  for (i in 1:length(comparison)) {
    object.net <- methods::slot(object, 'net')[[comparison[i]]]
    prob[[i]] = simplify2array(object.net)
  }
  names(prob) <- object.names
  
  if (is.null(point.shape)) {
    point.shape <- 1:15
  }
  
  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
    pathway.remove <- rownames(similarity)[which(colSums(similarity) == 1)]
    pathway.remove <- sub("--.*", "", pathway.remove)
  }
  
  if (length(pathway.remove) > 0) {
    for (i in 1:length(prob)) {
      probi <- prob[[i]]
      pathway.remove.idx <- which(paste0(dimnames(probi)[[3]],"--",object.names[i]) %in% pathway.remove)
      if (length(pathway.remove.idx) > 0) {
        probi <- probi[ , , -pathway.remove.idx]
      }
      prob[[i]] <- probi
    }
  }
  
  prob_sum.each <- list()
  signalingAll <- c()
  for (i in 1:length(prob)) {
    probi <- prob[[i]]
    prob_sum.each[[i]] <- apply(probi, 3, sum)
    signalingAll <- c(signalingAll, paste0(names(prob_sum.each[[i]]),"--",object.names[i]))
  }
  prob_sum <- unlist(prob_sum.each)
  names(prob_sum) <- signalingAll
  prob_sum <- prob_sum[rownames(Y)]
  
  group <- sub(".*--", "", names(prob_sum))
  labels = sub("--.*", "", names(prob_sum))
  
  df <- data.frame(x = Y[,1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum),
                   labels = as.character(labels), clusters = as.factor(clusters), group = factor(group, levels = unique(group)))
  
  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(clusters)))
  }
  
  ggAll <- vector("list", length(unique(clusters)))
  for (i in 1:length(unique(clusters))) {
    clusterID = i
    title <- paste0("Cluster ",  clusterID)
    df2 <- df[df$clusters %in% clusterID,]
    gg <- ggplot(data = df2, aes(x, y)) +
      geom_point(aes(size = Commun.Prob., shape = group), 
                 fill = alpha(color.use[clusterID], alpha = dot.alpha),
                 colour = alpha(color.use[clusterID], alpha = 1),
                 stroke = stroke.size) +
      CellChat_theme_opts() +
      theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in")) +
      guides(colour = guide_legend(override.aes = list(size = 3))) +
      labs(title = title, x = xlabel, y = ylabel) +
      scale_size_continuous(limits = c(0,1), range = dot.size, breaks = c(0.1,0.5,0.9)) +
      theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
      theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25)) +
      theme(plot.title = element_text(size = 50))
    idx <- match(unique(df2$group), levels(df$group), nomatch = 0)
    gg <- gg + scale_shape_manual(values = point.shape[idx])
    if (do.label) {
      gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels), 
                                          colour = color.use[clusterID], 
                                          size = label.size, 
                                          show.legend = F, 
                                          segment.size = 0.2, 
                                          segment.alpha = 0.5) +
        scale_alpha_discrete(range = c(1, 0.6))
    }
    
    if (!show.legend) {
      gg <- gg + theme(legend.position = "none")
    }
    
    if (!show.axes) {
      gg <- gg + theme_void()
    }
    ggAll[[i]] <- gg
  }
  gg.combined <- cowplot::plot_grid(plotlist = ggAll, ncol = nCol)
  
  gg.combined
}

#embeddingPairwiseZoomIn
p4 <- netVisual_embeddingPairwiseZoomIn_Neuron(neuronchat_merge_list, slot.name = "net_analysis", type = "functional", label.size = 20,comparison=c(1,2,3,4),nCol=3, dot.size = c(15, 20), stroke.size = 3)

ggsave(
  filename = "E:/master_graduation/neuronchatplot/embeddingPairwiseZoomIn.pdf",
  plot = p4,      
  device = "pdf",              
  width = 3840 * 3,           
  height = 2560 * 3,         
  units = "px",                 
  dpi = 300,
  limitsize = FALSE
)


#Heatmap
net1234 <- neuronchat_merge_list@net[c(1,2,3,4)]
net1 <- net1234[[1]];names(net1) <- paste(names(net1),'--td+ey-',sep='')
net2 <- net1234[[2]];names(net2) <- paste(names(net2),'--td-ey+',sep='')
net3 <- net1234[[3]];names(net3) <- paste(names(net3),'--td+ey+',sep='')
net4 <- net1234[[4]];names(net4) <- paste(names(net4),'--td-ey-',sep='')
net1234_list <- c(net1,net2,net3,net4)
interaction_group <- neuronchat_merge_list@net_analysis$similarity$functional$group$`1-2-3-4`
hlist <- list()
gb_heatmap <- list()
grid.newpage() 
x_seq <- c(0,0.2,0.4,0.6,0.8)
file_name <- map_output_path("E:/master_graduation/neuronchatplot/Heatmap.png")
dir.create(dirname(file_name), recursive = TRUE, showWarnings = FALSE)
png(filename = file_name, width = 17000, height = 8640, res = 300)  
for(j in 1:length(sort(unique(interaction_group),decreasing = F))){
  net_aggregated_group2 <- net_aggregation(net1234_list[names(interaction_group[interaction_group==j])],method = 'weight')
  col_map = brewer.pal(8,"YlOrBr");
  h <- Heatmap(net_aggregated_group2, 
               col = col_map,
               cluster_rows = FALSE,cluster_columns=FALSE,
               row_names_side='left',column_names_side='bottom',
               row_title='Sender',row_title_side='left',
               row_title_gp = gpar(fontsize = 40),
               column_title='Receiver',
               column_title_side = "bottom",
               column_title_gp = gpar(fontsize = 40),
               column_names_rot = 60,
               row_names_gp = gpar(fontsize = 40),
               column_names_gp = gpar(fontsize = 40),
               heatmap_legend_param = list(
                 title = "",
                 title_gp = gpar(fontsize = 40),
                 labels_gp = gpar(fontsize = 40),
                 grid_width = unit(1, "cm")
               ))
  gb_heatmap[[j]] = grid.grabExpr(draw(h,column_title=paste('pattern cluster',j), padding = unit(c(2, 2, 2, 2), "mm")))
  pushViewport(viewport(x = x_seq[j], y = 1, width = 0.2, height = 0.35, just = c("left", "top"),xscale = c(0, 1), yscale = c(0, 1)));grid.draw(gb_heatmap[[j]]);popViewport()
}
dev.off()
