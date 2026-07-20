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
# Differential expression and functional enrichment for all neuronal cells
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
library(cowplot)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(reshape2)
library(patchwork)
library(openxlsx)
library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db) 
library(enrichplot)
library(ReactomePA)



neurons <- readRDS(snakemake@input[["neurons"]])
neurons_con <- neurons %>%
  subset(group %in% grep("_con$", group, value = TRUE))

Idents(neurons_con) <- neurons_con@meta.data$group
acquisition_vs_other <- FindMarkers(neurons_con, ident.1 = "Acquisition_con", ident.2 = "Other_con")
retrieval_vs_other <- FindMarkers(neurons_con, ident.1 = "Retrieval_con", ident.2 = "Other_con")
overlapping_vs_other <- FindMarkers(neurons_con, ident.1 = "Overlapping_con", ident.2 = "Other_con")
overlapping_vs_acquisition <- FindMarkers(neurons_con, ident.1 = "Overlapping_con", ident.2 = "Acquisition_con")
overlapping_vs_retrieval <- FindMarkers(neurons_con, ident.1 = "Overlapping_con", ident.2 = "Retrieval_con")

deg_output_dir <- file.path(results_dir, "differential_expression")
dir.create(deg_output_dir, recursive = TRUE, showWarnings = FALSE)
deg_tables <- list(
  acquisition_vs_other = acquisition_vs_other,
  retrieval_vs_other = retrieval_vs_other,
  overlapping_vs_other = overlapping_vs_other,
  overlapping_vs_acquisition = overlapping_vs_acquisition,
  overlapping_vs_retrieval = overlapping_vs_retrieval
)
for (comparison_name in names(deg_tables)) {
  write.table(
    deg_tables[[comparison_name]],
    file = file.path(deg_output_dir, paste0(comparison_name, ".tsv")),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}


filter_DEGs <- function(deg_list) {
  deg_list_filtered <- deg_list %>%
    rownames_to_column(var = "gene") %>%  
    filter(abs(avg_log2FC) > 0.25 & p_val_adj < 0.05) %>%
    filter(!gene %in% c("tdTomato", "EYFP")) %>%  
    column_to_rownames(var = "gene")  
  return(deg_list_filtered)
}


acquisition_filtered <- filter_DEGs(acquisition_vs_other)
retrieval_filtered <- filter_DEGs(retrieval_vs_other)
overlapping_filtered <- filter_DEGs(overlapping_vs_other)
overlapping_vs_acquisition_filtered <- filter_DEGs(overlapping_vs_acquisition)
overlapping_vs_retrieval_filtered <- filter_DEGs(overlapping_vs_retrieval)


acquisition_genes <- rownames(acquisition_filtered)
retrieval_genes <- rownames(retrieval_filtered)
overlapping_genes <- rownames(overlapping_filtered)
overlapping_vs_acquisition_genes <- rownames(overlapping_vs_acquisition_filtered)
overlapping_vs_retrieval_genes <- rownames(overlapping_vs_retrieval_filtered)



run_GO_enrichment <- function(gene_list, ont_type = "BP") {
  go_results <- enrichGO(gene = gene_list,
                         OrgDb = org.Mm.eg.db,
                         keyType = "SYMBOL",
                         ont = ont_type,
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)
  return(go_results)
}


acquisition_go_bp <- run_GO_enrichment(acquisition_genes, "BP")
acquisition_go_mf <- run_GO_enrichment(acquisition_genes, "MF")
acquisition_go_cc <- run_GO_enrichment(acquisition_genes, "CC")

retrieval_go_bp <- run_GO_enrichment(retrieval_genes, "BP")
retrieval_go_mf <- run_GO_enrichment(retrieval_genes, "MF")
retrieval_go_cc <- run_GO_enrichment(retrieval_genes, "CC")

overlapping_go_bp <- run_GO_enrichment(overlapping_genes, "BP")
overlapping_go_mf <- run_GO_enrichment(overlapping_genes, "MF")
overlapping_go_cc <- run_GO_enrichment(overlapping_genes, "CC")

overlapping_vs_acquisition_go_bp <- run_GO_enrichment(overlapping_vs_acquisition_genes, "BP")
overlapping_vs_acquisition_go_mf <- run_GO_enrichment(overlapping_vs_acquisition_genes, "MF")
overlapping_vs_acquisition_go_cc <- run_GO_enrichment(overlapping_vs_acquisition_genes, "CC")

overlapping_vs_retrieval_go_bp <- run_GO_enrichment(overlapping_vs_retrieval_genes, "BP")
overlapping_vs_retrieval_go_mf <- run_GO_enrichment(overlapping_vs_retrieval_genes, "MF")
overlapping_vs_retrieval_go_cc <- run_GO_enrichment(overlapping_vs_retrieval_genes, "CC")



run_KEGG_enrichment <- function(gene_list) {
  gene_symbols <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  kegg_results <- enrichKEGG(gene = gene_symbols$ENTREZID,
                             organism = "mmu", 
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
  if (!is.null(kegg_results)) {
    kegg_results <- setReadable(kegg_results, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  }
  return(kegg_results)
}



acquisition_kegg <- run_KEGG_enrichment(acquisition_genes)
retrieval_kegg <- run_KEGG_enrichment(retrieval_genes)
overlapping_kegg <- run_KEGG_enrichment(overlapping_genes)
overlapping_vs_acquisition_kegg <- run_KEGG_enrichment(overlapping_vs_acquisition_genes)
overlapping_vs_retrieval_kegg <- run_KEGG_enrichment(overlapping_vs_retrieval_genes)



run_Reactome_enrichment <- function(gene_list) {
  gene_symbols <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  reactome_results <- enrichPathway(gene = gene_symbols$ENTREZID,
                                    organism = "mouse",
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.05)
  if (!is.null(reactome_results)) {
    reactome_results <- setReadable(reactome_results, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  }
  return(reactome_results)
}


acquisition_reactome <- run_Reactome_enrichment(acquisition_genes)
retrieval_reactome <- run_Reactome_enrichment(retrieval_genes)
overlapping_reactome <- run_Reactome_enrichment(overlapping_genes)
overlapping_vs_acquisition_reactome <- run_Reactome_enrichment(overlapping_vs_acquisition_genes)
overlapping_vs_retrieval_reactome <- run_Reactome_enrichment(overlapping_vs_retrieval_genes)


enrichment_list <- list(
  acquisition_go_bp = acquisition_go_bp,
  acquisition_go_mf = acquisition_go_mf,
  acquisition_go_cc = acquisition_go_cc,
  retrieval_go_bp = retrieval_go_bp,
  retrieval_go_mf = retrieval_go_mf,
  retrieval_go_cc = retrieval_go_cc,
  overlapping_go_bp = overlapping_go_bp,
  overlapping_go_mf = overlapping_go_mf,
  overlapping_go_cc = overlapping_go_cc,
  over_vs_acqui_go_bp = overlapping_vs_acquisition_go_bp,
  over_vs_acqui_go_mf = overlapping_vs_acquisition_go_mf,
  over_vs_acqui_go_cc = overlapping_vs_acquisition_go_cc,
  over_vs_retri_go_bp = overlapping_vs_retrieval_go_bp,
  over_vs_retri_go_mf = overlapping_vs_retrieval_go_mf,
  over_vs_retri_go_cc = overlapping_vs_retrieval_go_cc,
  acquisition_kegg = acquisition_kegg,
  retrieval_kegg = retrieval_kegg,
  overlapping_kegg = overlapping_kegg,
  over_vs_acqui_kegg = overlapping_vs_acquisition_kegg,
  over_vs_retri_kegg = overlapping_vs_retrieval_kegg,
  acquisition_reactome = acquisition_reactome,
  retrieval_reactome = retrieval_reactome,
  overlapping_reactome = overlapping_reactome,
  over_vs_acqui_reactome = overlapping_vs_acquisition_reactome,
  over_vs_retri_reactome = overlapping_vs_retrieval_reactome
)

wb <- createWorkbook()
for (name in names(enrichment_list)) {
  enrichment_result <- enrichment_list[[name]]
  if (!is.null(enrichment_result) && !is.null(enrichment_result@result)) {
    result_df <- enrichment_result@result
    addWorksheet(wb, sheetName = name)
    writeData(wb, sheet = name, x = result_df)
  } else {
    message("No result found for ", name)
  }
}
saveWorkbook(wb, file = "enrichment_results.xlsx", overwrite = TRUE)

# Selected functional-enrichment tables
library(openxlsx)  
library(dplyr)  


entry_list <- list(  
  acquisition_go_bp = c("GO:0031345", "GO:0010977", "GO:0050767"),  
  acquisition_go_mf = c("GO:0003779", "GO:0001227"),  
  acquisition_go_cc = c("GO:0045211", "GO:0097449"),  
  retrieval_go_bp = c("GO:0035235", "GO:1990806"),  
  retrieval_go_mf = c("GO:0019911"),  
  retrieval_go_cc = c("GO:0098802"),  
  overlapping_go_bp = c("GO:0001774", "GO:0045088"),  
  overlapping_go_mf = c("GO:0140375"),  
  overlapping_go_cc = c("GO:0062023"),  
  over_vs_acqui_go_bp = c("GO:0002449"),  
  over_vs_acqui_go_mf = c("GO:0140375"),  
  over_vs_acqui_go_cc = c("GO:0045335"),  
  over_vs_retri_go_bp = character(0),  
  over_vs_retri_go_mf = character(0),  
  over_vs_retri_go_cc = character(0),  
  acquisition_kegg = c("mmu04820", "mmu01200"),  
  retrieval_kegg = character(0),  
  overlapping_kegg = c("mmu05150"),  
  over_vs_acqui_kegg = c("mmu05150"),  
  over_vs_retri_kegg = character(0),  
  acquisition_reactome = c("R-MMU-1971475"),  
  retrieval_reactome = character(0),  
  overlapping_reactome = c("R-MMU-3000178"),  
  over_vs_acqui_reactome = c("R-MMU-977606"),  
  over_vs_retri_reactome = character(0)  
)  


enrichment_list <- list(  
  acquisition_go_bp = acquisition_go_bp,  
  acquisition_go_mf = acquisition_go_mf,  
  acquisition_go_cc = acquisition_go_cc,  
  retrieval_go_bp = retrieval_go_bp,  
  retrieval_go_mf = retrieval_go_mf,  
  retrieval_go_cc = retrieval_go_cc,  
  overlapping_go_bp = overlapping_go_bp,  
  overlapping_go_mf = overlapping_go_mf,  
  overlapping_go_cc = overlapping_go_cc,  
  over_vs_acqui_go_bp = overlapping_vs_acquisition_go_bp,  
  over_vs_acqui_go_mf = overlapping_vs_acquisition_go_mf,  
  over_vs_acqui_go_cc = overlapping_vs_acquisition_go_cc,  
  over_vs_retri_go_bp = overlapping_vs_retrieval_go_bp,  
  over_vs_retri_go_mf = overlapping_vs_retrieval_go_mf,  
  over_vs_retri_go_cc = overlapping_vs_retrieval_go_cc,  
  acquisition_kegg = acquisition_kegg,  
  retrieval_kegg = retrieval_kegg,  
  overlapping_kegg = overlapping_kegg,  
  over_vs_acqui_kegg = overlapping_vs_acquisition_kegg,  
  over_vs_retri_kegg = overlapping_vs_retrieval_kegg,  
  acquisition_reactome = acquisition_reactome,  
  retrieval_reactome = retrieval_reactome,  
  overlapping_reactome = overlapping_reactome,  
  over_vs_acqui_reactome = overlapping_vs_acquisition_reactome,  
  over_vs_retri_reactome = overlapping_vs_retrieval_reactome  
)  


wb <- createWorkbook()  


for (name in names(enrichment_list)) {  
  enrichment_result <- enrichment_list[[name]]  
  

  if (!is.null(enrichment_result) && !is.null(enrichment_result@result)) {  
    result_df <- as.data.frame(enrichment_result@result)
    

    if (name %in% names(entry_list)) {  
      term_ids <- entry_list[[name]]  
      filtered_data <- result_df[result_df$ID %in% term_ids, ]
      

      if (nrow(filtered_data) > 0) {  

        filtered_data <- filtered_data[, c("Description", "GeneRatio", "p.adjust", "geneID", "Count")]  
        colnames(filtered_data)[which(names(filtered_data) == "p.adjust")] <- "pvalue"  
        

        if ("GeneRatio" %in% colnames(filtered_data)) {  

          temp <- separate(filtered_data, GeneRatio, into = c("n", "N"), sep = "/", convert = TRUE)  

          temp$GeneRatio <- temp$n / temp$N  
          

          filtered_data <- temp[, c("Description", "GeneRatio", "pvalue", "geneID", "Count")]  
        }  
        
        addWorksheet(wb, sheetName = name)  
        writeData(wb, sheet = name, x = filtered_data)  
      } else {  
        message("No filtered result found for ", name)  
      }  
    } else {  
      message("No entry list found for ", name)  
    }  
  } else {  
    message("No result found for ", name)  
  }  
}  


saveWorkbook(wb, file = "enrichment_results_filtered.xlsx", overwrite = TRUE)


# Volcano plots for all neuronal comparisons
library(ggrepel)
library(ggplot2)
library(dplyr)
library(tibble)


volcano_plot <- function(DEG_list, objectname){
  DEG_list <- DEG_list %>%
    mutate(Difference = pct.1 - pct.2) %>%
    rownames_to_column("gene") %>%
    filter(!gene %in% c("tdTomato", "EYFP")) %>%
    mutate(threshold = case_when(
      avg_log2FC >= 0.25 & p_val_adj <= 0.05 ~ "Up",
      avg_log2FC <= -0.25 & p_val_adj <= 0.05 ~ "Down",
      TRUE ~ "false"
    ))

  plot <- ggplot(DEG_list, aes(x=Difference, y=avg_log2FC, color = threshold)) + 
          geom_point(size=3) + 
          scale_color_manual(values=c("Up" = "#C30078", "false" = "#CCCCCC","Down" = "#4370B4")) + 
          geom_text_repel(data=subset(DEG_list, avg_log2FC >= 0.25 & p_val_adj <= 0.05), 
                           aes(label=gene),
                           color="#C30078",
                           fill = NA,
                           label.padding = 0.1, 
                           size=7,
                           label.size = 0)+
          geom_text_repel(data=subset(DEG_list, avg_log2FC <= -0.25 & p_val_adj <= 0.05), 
                           aes(label=gene), 
                           color="#4370B4",
                           fill = NA,
                           label.padding = 0.1, 
                           size=7,
                           label.size = 0)+
          geom_vline(xintercept = 0.0,linetype=2)+
          geom_hline(yintercept = 0,linetype=2)+
          theme_classic()+
          theme(axis.title.x = element_text(size = 40),
                axis.title.y = element_text(size = 40),
                #legend.key.size = unit(2, "lines"),
                legend.text = element_text(size = 40),
                legend.title = element_text(size = 40),
                axis.text = element_text(size = 40))

  ggsave(
        filename = paste0("E:/master_graduation/volcano/pdf/", objectname, ".pdf"),
        plot = plot,      
        device = "pdf",              
        width = 3840,           
        height = 2160,         
        units = "px",                 
        dpi = 300,
        limitsize = FALSE
        )
}

DEG_list <- list(
  acquisition_vs_other = acquisition_vs_other,
  retrieval_vs_other = retrieval_vs_other,
  overlapping_vs_other = overlapping_vs_other,
  overlapping_vs_acquisition = overlapping_vs_acquisition,
  overlapping_vs_retrieval = overlapping_vs_retrieval
)

for (name in names(DEG_list)) {
  volcano_plot(DEG_list[[name]], name)
}

# Retrieval-versus-acquisition enrichment
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
library(cowplot)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(reshape2)
library(patchwork)
library(openxlsx)
library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db) 
library(enrichplot)
library(ReactomePA)



Idents(neurons_con) <- neurons_con@meta.data$group

retrieval_vs_acquisition <- FindMarkers(neurons_con, ident.1 = "Retrieval_con", ident.2 = "Acquisition_con")


filter_DEGs <- function(deg_list) {
  deg_list_filtered <- deg_list %>%
    rownames_to_column(var = "gene") %>%  
    filter(abs(avg_log2FC) > 0.25 & p_val_adj < 0.05) %>%
    filter(!gene %in% c("tdTomato", "EYFP")) %>%  
    column_to_rownames(var = "gene")  
  return(deg_list_filtered)
}



retrieval_vs_acquisition_filtered <- filter_DEGs(retrieval_vs_acquisition)


retrieval_vs_acquisition_genes <- rownames(retrieval_vs_acquisition_filtered)



run_GO_enrichment <- function(gene_list, ont_type = "BP") {
  go_results <- enrichGO(gene = gene_list,
                         OrgDb = org.Mm.eg.db,
                         keyType = "SYMBOL",
                         ont = ont_type,
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)
  return(go_results)
}


retrieval_vs_acquisition_go_bp <- run_GO_enrichment(retrieval_vs_acquisition_genes, "BP")
retrieval_vs_acquisition_go_mf <- run_GO_enrichment(retrieval_vs_acquisition_genes, "MF")
retrieval_vs_acquisition_go_cc <- run_GO_enrichment(retrieval_vs_acquisition_genes, "CC")



run_KEGG_enrichment <- function(gene_list) {
  gene_symbols <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  kegg_results <- enrichKEGG(gene = gene_symbols$ENTREZID,
                             organism = "mmu", 
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
  if (!is.null(kegg_results)) {
    kegg_results <- setReadable(kegg_results, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  }
  return(kegg_results)
}



retrieval_vs_acquisition_kegg <- run_KEGG_enrichment(retrieval_vs_acquisition_genes)



run_Reactome_enrichment <- function(gene_list) {
  gene_symbols <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  reactome_results <- enrichPathway(gene = gene_symbols$ENTREZID,
                                    organism = "mouse",
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.05)
  if (!is.null(reactome_results)) {
    reactome_results <- setReadable(reactome_results, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  }
  return(reactome_results)
}


retrieval_vs_acquisition_reactome <- run_Reactome_enrichment(retrieval_vs_acquisition_genes)


enrichment_list <- list(
  go_bp = retrieval_vs_acquisition_go_bp,
  go_mf = retrieval_vs_acquisition_go_mf,
  go_cc = retrieval_vs_acquisition_go_cc,
  kegg = retrieval_vs_acquisition_kegg,
  reactome = retrieval_vs_acquisition_reactome
)

wb <- createWorkbook()
for (name in names(enrichment_list)) {
  enrichment_result <- enrichment_list[[name]]
  if (!is.null(enrichment_result) && !is.null(enrichment_result@result)) {
    result_df <- enrichment_result@result
    addWorksheet(wb, sheetName = name)
    writeData(wb, sheet = name, x = result_df)
  } else {
    message("No result found for ", name)
  }
}
saveWorkbook(wb, file = "retrieval_vs_acquisition_enrichment_results.xlsx", overwrite = TRUE)


# Retrieval-versus-acquisition volcano plot
library(ggrepel)
library(ggplot2)
library(dplyr)
library(tibble)


volcano_plot <- function(DEG_list, objectname){
  DEG_list <- DEG_list %>%
    mutate(Difference = pct.1 - pct.2) %>%
    rownames_to_column("gene") %>%
    filter(!gene %in% c("tdTomato", "EYFP")) %>%
    mutate(threshold = case_when(
      avg_log2FC >= 0.25 & p_val_adj <= 0.05 ~ "Up",
      avg_log2FC <= -0.25 & p_val_adj <= 0.05 ~ "Down",
      TRUE ~ "false"
    ))
  
  plot <- ggplot(DEG_list, aes(x=Difference, y=avg_log2FC, color = threshold)) + 
    geom_point(size=3) + 
    scale_color_manual(values=c("Up" = "#C30078", "false" = "#CCCCCC","Down" = "#4370B4")) + 
    geom_text_repel(data=subset(DEG_list, avg_log2FC >= 0.25 & p_val_adj <= 0.05), 
                     aes(label=gene),
                     color="#C30078",
                     fill = NA,
                     label.padding = 0.1, 
                     size=7,
                     label.size = 0)+
    geom_text_repel(data=subset(DEG_list, avg_log2FC <= -0.25 & p_val_adj <= 0.05), 
                     aes(label=gene), 
                     color="#4370B4",
                     fill = NA,
                     label.padding = 0.1, 
                     size=7,
                     label.size = 0)+
    geom_vline(xintercept = 0.0,linetype=2)+
    geom_hline(yintercept = 0,linetype=2)+
    theme_classic()+
    theme(axis.title.x = element_text(size = 40),
          axis.title.y = element_text(size = 40),
          legend.key.size = unit(2, "lines"),
          legend.text = element_text(size = 40),
          legend.title = element_text(size = 40),
          axis.text = element_text(size = 40))
  
  ggsave(
    filename = paste0("E:/master_graduation/volcano/pdf/", objectname, ".pdf"),
    plot = plot,      
    device = "pdf",              
    width = 3840,           
    height = 2160,         
    units = "px",                 
    dpi = 300,
    limitsize = FALSE
  )
}

DEG_list <- list(
  retrieval_vs_acquisition = retrieval_vs_acquisition
)

for (name in names(DEG_list)) {
  volcano_plot(DEG_list[[name]], name)
}

# Retrieval-versus-acquisition enrichment plots
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
library(cowplot)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(reshape2)
library(patchwork)
library(openxlsx)
library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db) 
library(enrichplot)
library(ReactomePA)



extract_gene_pathway <- function(df) {
  gene_list <- strsplit(df$geneID, "/")
  result <- data.frame(
    Gene = unlist(gene_list),
    Pathway = rep(df$Description, lengths(gene_list))
  )
  result$Pathway <- gsub(" - Mus musculus \\(house mouse\\)", "", result$Pathway)
  return(result)
}

wrap_pathway <- function(text) {  
  words <- str_extract_all(text, "\\w+(?:[.-]\\w+)*")[[1]]
  if (length(words) <= 3) {  
    return(text)
  }  
  

  wrapped <- str_c(head(words, 3), collapse = " ")
  for (i in seq(from = 4, to = length(words), by = 3)) {  
    wrapped <- str_c(wrapped, "\n", str_c(words[i:min(i + 2, length(words))], collapse = " "))  
  }  
  return(wrapped)  
}  


custom_round_labels <- function(x) {  

  rounded_labels <- x  
  
  for (i in seq_along(x)) {  

    if (grepl("\\.", x[i])) {  
      decimal_part <- nchar(sub(".*\\.", "", as.character(x[i])))
      if (decimal_part >= 3) {  
        rounded_labels[i] <- round(x[i], 2)
      }  
    }  
  }  
  return(rounded_labels)  
}  

generate_enrichment_plots <- function(enrichment_result, filename_prefix) {
  if (is.null(enrichment_result) || nrow(enrichment_result@result) == 0) {
    message("No results for ", filename_prefix)
    return()
  }
  
  entry_ids <- entry_list[[filename_prefix]]
  
  if (is.null(entry_ids) || length(entry_ids) == 0) {
    message("No entry IDs defined for ", filename_prefix)
    return()
  } else {

    specific_entries <- enrichment_result@result %>%
      filter(ID %in% entry_ids)

    if (nrow(specific_entries) == 0) {
      message("No matching entries found for ", filename_prefix)
      return()
    }
  }
  

  specific_entries <- specific_entries %>%
    separate(GeneRatio, into = c("n", "N"), sep = "/", convert = TRUE) %>%
    mutate(GeneRatio = n / N)
  

  # top10 <- enrichment_result@result %>%
  # arrange(pvalue) %>%
  # head(10) %>%
  # separate(GeneRatio, into = c("n", "N"), sep = "/", convert = TRUE) %>%
  # mutate(GeneRatio = n / N)
  #if (nrow(top10) == 0) {
  # message("No terms in ", filename_prefix)
  # return()
  #}
  

  sankey_data <- extract_gene_pathway(specific_entries)
  sankey_data <- sankey_data %>%  
    group_by(Pathway) %>%  
    slice_head(n = 10) %>%
    ungroup()  
  sankey_data <- sankey_data %>%  
    mutate(Pathway = sapply(Pathway, wrap_pathway))
  
  specific_entries$Description <- factor(specific_entries$Description,levels = rev(specific_entries$Description))
  
  specific_entries_2 <- specific_entries[nrow(specific_entries):1, ]
  specific_entries_2 <- specific_entries_2 %>%
    mutate(ymax = cumsum(Count * 2)) %>%
    mutate(ymin = ymax -Count * 2) %>%
    mutate(label = (ymin + ymax)/2)
  
  mytheme<- theme(axis.title = element_text(size = 60),
                  axis.text = element_text(size = 60),
                  legend.title = element_text(size = 60, margin = margin(b = 20)),
                  legend.text = element_text(size = 60),
                  legend.key.height = unit(2.5, "cm"),
                  legend.margin = margin(b = 20),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank())
  
  p1<- ggplot() +
    geom_point(data = specific_entries_2,
               aes(x = -log10(pvalue),
                   y= label,
                   size= Count,
                   color= GeneRatio)) +
    scale_size_continuous(range=c(10,20)) +
    scale_colour_distiller(palette = "Reds", direction = 1) +
    labs(x = "-log10(pvalue)",
         y= "") +
    theme_bw() +
    mytheme +
    guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2))
  
  if (nrow(specific_entries) == 1) {  
    p1 <- p1 + scale_x_continuous(limits = c(-log10(specific_entries_2$pvalue)-0.5, -log10(specific_entries_2$pvalue)+0.5),   
                                  breaks = seq(-log10(specific_entries_2$pvalue)-0.5, -log10(specific_entries_2$pvalue)+0.5, by = 0.5),
                                  labels = function(x) format(round(x, 1), nsmall = 1))   
  } else {  
    p1 <- p1 + scale_x_continuous(labels = custom_round_labels)
  }  
  
  df_sankey <- sankey_data %>%
    make_long(Gene, Pathway)
  
  df_sankey$node <- factor(df_sankey$node,levels = c(sankey_data$Pathway %>% unique%>% rev,
                                                     sankey_data$Gene %>% unique %>% rev))
  
  mycol_up<- c4a('rainbow_wh_rd',length(unique(df_sankey$node)))
  
  p2<- ggplot(df_sankey, aes(x = x,
                             next_x= next_x,
                             node= node,
                             next_node= next_node,
                             fill= node,
                             label= node)) +
    geom_alluvial(flow.alpha = 0.5,
                  flow.fill = 'grey',
                  flow.color = 'grey80',
                  node.fill = mycol_up,
                  smooth= 8,
                  width= 0.05) +
    geom_alluvial_text(aes(hjust = 1),
                       size = 21,
                       position = position_nudge(x = -0.04),
                       color= "black")+
    theme_void() +
    theme(legend.position = 'none')
  
  p3 <- p2 + theme(plot.margin = unit(c(0,10,0,0),units= "cm"))
  
  plot <- ggdraw() + draw_plot(p3) + draw_plot(p1, scale = 0.67, x = 0.61, y=-0.23, width=0.43, height=1.39)
  
  ggsave(
    filename = paste0("E:/master_graduation/enrichplot/pdf/", filename_prefix, "_enrichment.pdf"),
    plot = plot,      
    device = "pdf",              
    width = 3840 * 3,           
    height = 2560 * 3,         
    units = "px",                 
    dpi = 300,
    limitsize = FALSE
  )
}


enrichment_list <- list(
  retrieval_vs_acquisition_go_bp = retrieval_vs_acquisition_go_bp,
  retrieval_vs_acquisition_go_mf = retrieval_vs_acquisition_go_mf,
  retrieval_vs_acquisition_go_cc = retrieval_vs_acquisition_go_cc,
  retrieval_vs_acquisition_kegg = retrieval_vs_acquisition_kegg,
  retrieval_vs_acquisition_reactome = retrieval_vs_acquisition_reactome
)

entry_list <- list(
  retrieval_vs_acquisition_go_bp = c("GO:0042552", "GO:0007272", "GO:0008366", "GO:0060249", "GO:0042063", "GO:0070444", "GO:0010959"),
  retrieval_vs_acquisition_go_mf = c("GO:0019911", "GO:0005506", "GO:0005516", "GO:0016783"),
  retrieval_vs_acquisition_go_cc = c("GO:0043209", "GO:0098802", "GO:0051286", "GO:0060187", "GO:0043218", "GO:0005796", "GO:0098636"),
  retrieval_vs_acquisition_kegg = c("mmu00790", "mmu04216", "mmu04978", "mmu04350", "mmu04066"),
  retrieval_vs_acquisition_reactome = c("R-MMU-196854", "R-MMU-2022928", "R-MMU-917937")
)

for (name in names(enrichment_list)) {
  generate_enrichment_plots(enrichment_list[[name]], name)
}

# Labeled enrichment plots
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
library(cowplot)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(reshape2)
library(patchwork)
library(openxlsx)
library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db) 
library(enrichplot)
library(ReactomePA)
library(scales)



extract_gene_pathway <- function(df) {
  gene_list <- strsplit(df$geneID, "/")
  result <- data.frame(
    Gene = unlist(gene_list),
    Pathway = rep(df$Description, lengths(gene_list))
  )
  result$Pathway <- gsub(" - Mus musculus \\(house mouse\\)", "", result$Pathway)
  #result$Pathway <- gsub(", RNA polymerase II-specific", "", result$Pathway)
  return(result)
}

wrap_pathway <- function(text) {  
  words <- str_extract_all(text, "\\w+(?:[.-]\\w+)*")[[1]]
  if (length(words) <= 3) {  
    return(text)
  }  
  

  wrapped <- str_c(head(words, 3), collapse = " ")
  for (i in seq(from = 4, to = length(words), by = 3)) {  
    wrapped <- str_c(wrapped, "\n", str_c(words[i:min(i + 2, length(words))], collapse = " "))  
  }  
  return(wrapped)  
}  


custom_round_labels <- function(x) {  

  rounded_labels <- x  
  
  for (i in seq_along(x)) {  

    if (grepl("\\.", x[i])) {  
      decimal_part <- nchar(sub(".*\\.", "", as.character(x[i])))
      if (decimal_part >= 3) {  
        rounded_labels[i] <- round(x[i], 2)
      }  
    }  
  }  
  return(rounded_labels)  
}  

generate_enrichment_plots <- function(enrichment_result, filename_prefix) {
  if (is.null(enrichment_result) || nrow(enrichment_result@result) == 0) {
    message("No results for ", filename_prefix)
    return()
  }
  
  entry_ids <- entry_list[[filename_prefix]]

  if (is.null(entry_ids) || length(entry_ids) == 0) {
    message("No entry IDs defined for ", filename_prefix)
    return()
  } else {

    specific_entries <- enrichment_result@result %>%
      filter(ID %in% entry_ids)

    if (nrow(specific_entries) == 0) {
      message("No matching entries found for ", filename_prefix)
      return()
    }
  }
  

  specific_entries <- specific_entries %>%
    separate(GeneRatio, into = c("n", "N"), sep = "/", convert = TRUE) %>%
    mutate(GeneRatio = n / N)
  

  # top10 <- enrichment_result@result %>%
  # arrange(pvalue) %>%
  # head(10) %>%
  # separate(GeneRatio, into = c("n", "N"), sep = "/", convert = TRUE) %>%
  # mutate(GeneRatio = n / N)
  #if (nrow(top10) == 0) {
  # message("No terms in ", filename_prefix)
  # return()
  #}
  

  sankey_data <- extract_gene_pathway(specific_entries)
  sankey_data <- sankey_data %>%  
    group_by(Pathway) %>%  
    slice_head(n = 10) %>%
    ungroup()  
  sankey_data <- sankey_data %>%  
    mutate(Pathway = sapply(Pathway, wrap_pathway))
  print(head(sankey_data))
  
  df_sankey <- sankey_data %>%
    make_long(Gene, Pathway)
  
  df_sankey$node <- factor(df_sankey$node,levels = c(sankey_data$Pathway %>% unique%>% rev,
                                                     sankey_data$Gene %>% unique %>% rev))
  print(head(df_sankey))
  
  specific_entries$Description <- factor(specific_entries$Description,levels = rev(specific_entries$Description))
  
  specific_entries_2 <- specific_entries[nrow(specific_entries):1, ]
  specific_entries_2 <- specific_entries_2 %>%
    mutate(ymax = cumsum(Count * 2)) %>%
    mutate(ymin = ymax -Count * 2) %>%
    mutate(label = (ymin + ymax)/2)
  
  mytheme<- theme(axis.title = element_text(size = 60),
                  axis.text = element_text(size = 60),
                  legend.title = element_text(size = 60, margin = margin(b = 20)),
                  legend.text = element_text(size = 60),
                  legend.key.height = unit(2.5, "cm"),
                  legend.margin = margin(b = 20),
                  plot.margin = unit(c(0, 0, 0, 0), "cm"),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank())
  
  p1<- ggplot() +
    geom_point(data = specific_entries_2,
               aes(x = -log10(pvalue),
                   y= label,
                   size= Count,
                   color= GeneRatio)) +
    scale_size_continuous(range=c(10,20)) +
    scale_colour_distiller(palette = "Reds", direction = 1, labels = function(x) format(x, nsmall = 3, digits = 3)) +
    #scale_colour_gradient(labels = function(x) format(x, nsmall = 4, digits = 4)) +  
    labs(x = expression(-log[10](p))) +
    theme_bw() +
    mytheme +
    guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2))
  
  if (nrow(specific_entries) == 1) {  
    p1 <- p1 + scale_x_continuous(limits = c(-log10(specific_entries_2$pvalue)-0.5, -log10(specific_entries_2$pvalue)+0.5),   
                                  breaks = seq(-log10(specific_entries_2$pvalue)-0.5, -log10(specific_entries_2$pvalue)+0.5, by = 0.5),
                                  labels = function(x) format(round(x, 1), nsmall = 1))   
  } else {  
    p1 <- p1 + scale_x_continuous(labels = custom_round_labels)
  }  
  
  
  
  mycol_up<- c4a('rainbow_wh_rd',length(unique(df_sankey$node)))
  
  p2<- ggplot(df_sankey, aes(x = x,
                                next_x= next_x,
                                node= node,
                                next_node= next_node,
                                fill= node,
                                label= node)) +
    geom_alluvial(flow.alpha = 0.5,
                flow.fill = 'grey',
                flow.color = 'grey80',
                node.fill = mycol_up,
                smooth= 8,
                width= 0.05) +
    geom_alluvial_text(aes(hjust = 1),
                     size = 21,
                     position = position_nudge(x = -0.04),
                     color= "black")+
    theme_void() +
    theme(legend.position = 'none')
  
  p3 <- p2 + theme(plot.margin = unit(c(0,10,0,0),units= "cm"))
  
  plot <- ggdraw() + draw_plot(p3) + draw_plot(p1, scale = 0.67, x = 0.61, y=-0.23, width=0.43, height=1.39)
  
  ggsave(
    filename = paste0("E:/master_graduation/enrichplot/pdf/", filename_prefix, "_enrichment.pdf"),
    plot = plot,      
    device = "pdf",              
    width = 3840 * 3,           
    height = 2560 * 3,         
    units = "px",                 
    dpi = 300,
    limitsize = FALSE
  )
}


enrichment_list <- list(
  acquisition_go_bp = acquisition_go_bp,
  acquisition_go_mf = acquisition_go_mf,
  acquisition_go_cc = acquisition_go_cc,
  retrieval_go_bp = retrieval_go_bp,
  retrieval_go_mf = retrieval_go_mf,
  retrieval_go_cc = retrieval_go_cc,
  overlapping_go_bp = overlapping_go_bp,
  overlapping_go_mf = overlapping_go_mf,
  overlapping_go_cc = overlapping_go_cc,
  overlapping_vs_acquisition_go_bp = overlapping_vs_acquisition_go_bp,
  overlapping_vs_acquisition_go_mf = overlapping_vs_acquisition_go_mf,
  overlapping_vs_acquisition_go_cc = overlapping_vs_acquisition_go_cc,
  overlapping_vs_retrieval_go_bp = overlapping_vs_retrieval_go_bp,
  overlapping_vs_retrieval_go_mf = overlapping_vs_retrieval_go_mf,
  overlapping_vs_retrieval_go_cc = overlapping_vs_retrieval_go_cc,
  acquisition_kegg = acquisition_kegg,
  retrieval_kegg = retrieval_kegg,
  overlapping_kegg = overlapping_kegg,
  overlapping_vs_acquisition_kegg = overlapping_vs_acquisition_kegg,
  overlapping_vs_retrieval_kegg = overlapping_vs_retrieval_kegg,
  acquisition_reactome = acquisition_reactome,
  retrieval_reactome = retrieval_reactome,
  overlapping_reactome = overlapping_reactome,
  overlapping_vs_acquisition_reactome = overlapping_vs_acquisition_reactome,
  overlapping_vs_retrieval_reactome = overlapping_vs_retrieval_reactome
)

entry_list <- list(
  acquisition_go_bp = c("GO:0031345", "GO:0010977", "GO:0050767"),
  acquisition_go_mf = c("GO:0003779", "GO:0001227"),
  acquisition_go_cc = c("GO:0045211", "GO:0097449", "GO:0097386"),
  retrieval_go_bp = c("GO:0035235", "GO:1990806", "GO:0007272", "GO:0008366", "GO:0007215"),
  retrieval_go_mf = c("GO:0019911", "GO:0005516"),
  retrieval_go_cc = c("GO:0098802"),
  overlapping_go_bp = c("GO:0001774", "GO:0045088", "GO:0002269", "GO:0032760"),
  overlapping_go_mf = c("GO:0140375"),
  overlapping_go_cc = c("GO:0062023", "GO:0043235", "GO:0005581"),
  overlapping_vs_acquisition_go_bp = c("GO:0002449", "GO:0001906", "GO:0098883", "GO:0001909"),
  overlapping_vs_acquisition_go_mf = c("GO:0140375", "GO:0016505", "GO:0023026"),
  overlapping_vs_acquisition_go_cc = c("GO:0045335", "GO:0042613"),
  overlapping_vs_retrieval_go_bp = character(0),  
  overlapping_vs_retrieval_go_mf = character(0),
  overlapping_vs_retrieval_go_cc = character(0),
  acquisition_kegg = c("mmu01200", "mmu01230", "mmu04922"),
  retrieval_kegg = character(0),
  overlapping_kegg = c("mmu05152", "mmu05140","mmu04145", "mmu05323"),
  overlapping_vs_acquisition_kegg = c("mmu05150", "mmu05140", "mmu05322", "mmu04145", "mmu05171"),
  overlapping_vs_retrieval_kegg = character(0),
  acquisition_reactome = c("R-MMU-71387", "R-MMU-77289", "R-MMU-1630316"),
  retrieval_reactome = character(0),
  overlapping_reactome = c("R-MMU-3000178", "R-MMU-5686938", "R-MMU-977606", "R-MMU-1474244", "R-MMU-166658", "R-MMU-166663"),
  overlapping_vs_acquisition_reactome = c("R-MMU-977606", "R-MMU-166663", "R-MMU-166658", "R-MMU-166786", "R-MMU-1222556", "R-MMU-512988", "R-MMU-2029481"),
  overlapping_vs_retrieval_reactome = character(0)
)

for (name in names(enrichment_list)) {
  generate_enrichment_plots(enrichment_list[[name]], name)
}

# Unlabeled enrichment plots
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
library(cowplot)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(reshape2)
library(patchwork)
library(openxlsx)
library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db) 
library(enrichplot)
library(ReactomePA)



extract_gene_pathway <- function(df) {
  gene_list <- strsplit(df$geneID, "/")
  result <- data.frame(
    Gene = unlist(gene_list),
    Pathway = rep(df$Description, lengths(gene_list))
  )
  result$Pathway <- gsub(" - Mus musculus \\(house mouse\\)", "", result$Pathway)
  result$Pathway <- gsub(", RNA polymerase II-specific", "", result$Pathway)
  return(result)
}

generate_enrichment_plots <- function(enrichment_result, filename_prefix) {
  if (is.null(enrichment_result) || nrow(enrichment_result@result) == 0) {
    message("No results for ", filename_prefix)
    return()
  }
  
  entry_ids <- entry_list[[filename_prefix]]
  
  if (is.null(entry_ids) || length(entry_ids) == 0) {
    message("No entry IDs defined for ", filename_prefix)
    return()
  } else {

    specific_entries <- enrichment_result@result %>%
      filter(ID %in% entry_ids)

    if (nrow(specific_entries) == 0) {
      message("No matching entries found for ", filename_prefix)
      return()
    }
  }
  

  specific_entries <- specific_entries %>%
    separate(GeneRatio, into = c("n", "N"), sep = "/", convert = TRUE) %>%
    mutate(GeneRatio = n / N)
  

  # top10 <- enrichment_result@result %>%
  # arrange(pvalue) %>%
  # head(10) %>%
  # separate(GeneRatio, into = c("n", "N"), sep = "/", convert = TRUE) %>%
  # mutate(GeneRatio = n / N)
  #if (nrow(top10) == 0) {
  # message("No terms in ", filename_prefix)
  # return()
  #}
  

  sankey_data <- extract_gene_pathway(specific_entries)
  print(head(sankey_data))
  
  df_sankey <- sankey_data %>%
    make_long(Gene, Pathway)
  
  df_sankey$node <- factor(df_sankey$node,levels = c(sankey_data$Pathway %>% unique%>% rev,
                                                     sankey_data$Gene %>% unique %>% rev))
  print(head(df_sankey))
  
  specific_entries$Description <- factor(specific_entries$Description,levels = rev(specific_entries$Description))
  
  specific_entries_2 <- specific_entries[nrow(specific_entries):1, ]
  specific_entries_2 <- specific_entries_2 %>%
    mutate(ymax = cumsum(Count * 2)) %>%
    mutate(ymin = ymax -Count * 2) %>%
    mutate(label = (ymin + ymax)/2)
  
  mytheme<- theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  legend.title = element_blank(),
                  legend.text = element_blank(),
                  legend.key.height = unit(2.5, "cm"),
                  legend.margin = margin(b = 40),
                  plot.margin = unit(c(0, 0, 0, 0), "cm"),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank())
  
  p1<- ggplot() +
    geom_point(data = specific_entries_2,
               aes(x = -log10(pvalue),
                   y= label,
                   size= Count,
                   color= GeneRatio)) +
    scale_size_continuous(range=c(10,20)) +
    scale_colour_distiller(palette = "Reds", direction = 1) +
    labs(x = "-log10(pvalue)",
         y= "") +
    theme_bw() +
    mytheme
  
  
  
  mycol_up<- c4a('rainbow_wh_rd',length(unique(df_sankey$node)))
  
  p2<- ggplot(df_sankey, aes(x = x,
                             next_x= next_x,
                             node= node,
                             next_node= next_node,
                             fill= node,
                             label= node)) +
    geom_alluvial(flow.alpha = 0.5,
                  flow.fill = 'grey',
                  flow.color = 'grey80',
                  node.fill = mycol_up,
                  smooth= 8,
                  width= 0.08) +

                       #size = 10,
                       #position = position_nudge(x = -0.05),
                       #color= "black")+
    theme_void() +
    theme(legend.position = 'none')
  
  p3 <- p2 + theme(plot.margin = unit(c(0,10,0,0),units= "cm"))
  
  plot <- ggdraw() + draw_plot(p3) + draw_plot(p1, scale = 0.67, x = 0.61, y=-0.21, width=0.43, height=1.39)
  
  ggsave(
    filename = paste0("E:/master_graduation/enrichplot/nolabel/", filename_prefix, "_enrichment_nolabel.png"),
    plot = plot,      
    device = "png",              
    width = 3840 * 3,           
    height = 2560 * 4,         
    units = "px",                 
    dpi = 300,
    limitsize = FALSE
  )
}


enrichment_list <- list(
  acquisition_go_bp = acquisition_go_bp,
  acquisition_go_mf = acquisition_go_mf,
  acquisition_go_cc = acquisition_go_cc,
  retrieval_go_bp = retrieval_go_bp,
  retrieval_go_mf = retrieval_go_mf,
  retrieval_go_cc = retrieval_go_cc,
  overlapping_go_bp = overlapping_go_bp,
  overlapping_go_mf = overlapping_go_mf,
  overlapping_go_cc = overlapping_go_cc,
  overlapping_vs_acquisition_go_bp = overlapping_vs_acquisition_go_bp,
  overlapping_vs_acquisition_go_mf = overlapping_vs_acquisition_go_mf,
  overlapping_vs_acquisition_go_cc = overlapping_vs_acquisition_go_cc,
  overlapping_vs_retrieval_go_bp = overlapping_vs_retrieval_go_bp,
  overlapping_vs_retrieval_go_mf = overlapping_vs_retrieval_go_mf,
  overlapping_vs_retrieval_go_cc = overlapping_vs_retrieval_go_cc,
  acquisition_kegg = acquisition_kegg,
  retrieval_kegg = retrieval_kegg,
  overlapping_kegg = overlapping_kegg,
  overlapping_vs_acquisition_kegg = overlapping_vs_acquisition_kegg,
  overlapping_vs_retrieval_kegg = overlapping_vs_retrieval_kegg,
  acquisition_reactome = acquisition_reactome,
  retrieval_reactome = retrieval_reactome,
  overlapping_reactome = overlapping_reactome,
  overlapping_vs_acquisition_reactome = overlapping_vs_acquisition_reactome,
  overlapping_vs_retrieval_reactome = overlapping_vs_retrieval_reactome
)

entry_list <- list(
  acquisition_go_bp = c("GO:0031345", "GO:0010977"),
  acquisition_go_mf = c("GO:0003779", "GO:0001227"),
  acquisition_go_cc = c("GO:0045211", "GO:0097449", "GO:0097386"),
  retrieval_go_bp = c("GO:0035235", "GO:1990806", "GO:0007272", "GO:0008366", "GO:0007215"),
  retrieval_go_mf = c("GO:0019911", "GO:0005516"),
  retrieval_go_cc = c("GO:0098802"),
  overlapping_go_bp = c("GO:0001774", "GO:0045088", "GO:0002269", "GO:0032760"),
  overlapping_go_mf = c("GO:0140375"),
  overlapping_go_cc = c("GO:0062023", "GO:0043235", "GO:0005581"),
  overlapping_vs_acquisition_go_bp = c("GO:0002449", "GO:0001906", "GO:0098883", "GO:0001909"),
  overlapping_vs_acquisition_go_mf = c("GO:0140375", "GO:0016505", "GO:0023026"),
  overlapping_vs_acquisition_go_cc = c("GO:0045335", "GO:0042613"),
  overlapping_vs_retrieval_go_bp = character(0),  
  overlapping_vs_retrieval_go_mf = character(0),
  overlapping_vs_retrieval_go_cc = character(0),
  acquisition_kegg = c("mmu01200", "mmu01230", "mmu04922"),
  retrieval_kegg = character(0),
  overlapping_kegg = c("mmu05150", "mmu05152", "mmu05140", "mmu04145", "mmu05323"),
  overlapping_vs_acquisition_kegg = c("mmu05150", "mmu05140", "mmu05322", "mmu04145", "mmu05171"),
  overlapping_vs_retrieval_kegg = character(0),
  acquisition_reactome = c("R-MMU-71387", "R-MMU-77289", "R-MMU-1630316"),
  retrieval_reactome = character(0),
  overlapping_reactome = c("R-MMU-3000178", "R-MMU-5686938", "R-MMU-977606", "R-MMU-1474244", "R-MMU-500792", "R-MMU-166658", "R-MMU-166663"),
  overlapping_vs_acquisition_reactome = c("R-MMU-977606", "R-MMU-166663", "R-MMU-166658", "R-MMU-166786", "R-MMU-1222556", "R-MMU-512988", "R-MMU-2029481"),
  overlapping_vs_retrieval_reactome = character(0)
)

for (name in names(enrichment_list)) {
  generate_enrichment_plots(enrichment_list[[name]], name)
}

# Marker-expression heatmap
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(grid)
library(Seurat)
library(reshape2)

Idents(neurons_con) <- neurons_con@meta.data$group
Acquisition_character <- FindMarkers(neurons_con, ident.1 = "Acquisition_con", logfc.threshold = 0.25)
Retrieval_character <- FindMarkers(neurons_con, ident.1 = "Retrieval_con", logfc.threshold = 0.25)
Overlapping_character <- FindMarkers(neurons_con, ident.1 = "Overlapping_con", logfc.threshold = 0.25)

filter_genes <- function(df) {
  df_filtered <- df %>%
    filter(abs(avg_log2FC) > 1 & p_val_adj < 0.05)
  return(df_filtered)
}

Acquisition_character_filtered <- filter_genes(Acquisition_character)
Retrieval_character_filtered <- filter_genes(Retrieval_character)
Overlapping_character_filtered <- filter_genes(Overlapping_character)

get_top10_by_log2FC <- function(df) {
  df %>%
    arrange(desc(abs(avg_log2FC))) %>% 
    head(10)
}

top10_Acquisition_character <- get_top10_by_log2FC(Acquisition_character_filtered)
top10_Retrieval_character <- get_top10_by_log2FC(Retrieval_character_filtered)
top10_Overlapping_character <- get_top10_by_log2FC(Overlapping_character_filtered)

combined_top10_character <- unique(c(rownames(top10_Acquisition_character), 
                                 rownames(top10_Retrieval_character), 
                                 rownames(top10_Overlapping_character)))

expression_matrix <- GetAssayData(neurons_con, layer = "data")
meta_data <- neurons_con@meta.data
calculate_group_mean <- function(genes, meta_data, group_name) {
  cells_in_group <- rownames(meta_data[meta_data$group == group_name, ])
  sub_matrix <- expression_matrix[genes, cells_in_group, drop = FALSE]
  rowMeans(as.matrix(sub_matrix))
}

character_means <- data.frame(
  Other_con = calculate_group_mean(combined_top10_character, meta_data, "Other_con"),
  Acquisition_con = calculate_group_mean(combined_top10_character, meta_data, "Acquisition_con"),
  Retrieval_con = calculate_group_mean(combined_top10_character, meta_data, "Retrieval_con"),
  Overlapping_con = calculate_group_mean(combined_top10_character, meta_data, "Overlapping_con")
)

character_matrix <- as.matrix(character_means)
rownames(character_matrix) <- combined_top10_character
colnames(character_matrix) <- c("td-ey-", "td+ey-", "td-ey+", "td+ey+")

colors <- c("#FFFFFF", "yellow", "red")

convert_to_long <- function(matrix, group_name) {
  df <- as.data.frame(matrix)
  df$Gene <- rownames(matrix)
  df_long <- melt(df, id.vars = "Gene", variable.name = "Group", value.name = "Expression")
  df_long$Group <- factor(df_long$Group, levels = colnames(matrix))
  return(df_long)
}

character_matrix_long <- convert_to_long(character_matrix)

plot_heatmap <- function(data, colors) {
  data$Gene <- factor(data$Gene, levels = (unique(data$Gene)))
  

  x_min <- 0.5
  x_max <- length(unique(data$Group)) + 0.5
  y_min <- 0.5
  y_max <- length(unique(data$Gene)) + 0.5
  
  ggplot(data, aes(x = Group, y = Gene)) +
    geom_point(aes(size = abs(Expression), fill = Expression), shape = 21, stroke = 0.5) + 
    geom_rect(
      xmin = x_min, xmax = x_max,
      ymin = y_min, ymax = y_max,
      color = "black", fill = NA, linewidth = 1
    ) +
    scale_fill_gradient2(
      low = colors[1], mid = colors[2], high = colors[3], midpoint = 0,
      name = "Expression"
    ) +
    scale_size_continuous(range = c(5, 15), name = "Expression") +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 40, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 40, hjust = 1),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.text = element_text(size = 40),
      legend.title = element_text(size = 40, margin = margin(b = 20)),
      legend.margin = margin(b = 20),
      legend.key.height = unit(2, "cm")
      #legend.key.size = unit(4, "lines")
    ) +
    coord_fixed(ratio = 1)
}

p <- plot_heatmap(character_matrix_long, colors)

ggsave(
  filename = "E:/master_graduation/heatmap/pdf/heatmap.pdf",
  plot = p,      
  device = "pdf",              
  width = 3840 * 2,           
  height = 2560 * 2,         
  units = "px",                 
  dpi = 300,
  limitsize = FALSE
)
