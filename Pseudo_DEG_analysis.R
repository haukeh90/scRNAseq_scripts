library(Libra)
library(dplyr)
library(Seurat)
library(openxlsx)
library(readr)

#Read in dataset----------------------------------------------------------------
seurat_object <- readRDS("~/directory/seurat_object.RDS")
table(seurat_object$dataset, seurat_object$bio_replicate)

#set default assay to RNA3 
DefaultAssay(seurat_object) <- "RNA3"
Idents(seurat_object) <- "final_clusters"
levels(seurat_object)

#Perform pseudo bulk analysis with LIRBA----------------------------------------

#check on relevant metadata 
unique(seurat_object$final_clusters)
unique(seurat_object$bio_replicate)
unique(seurat_object$Group)

#set correct factor for FC - check afterwards 
Idents(seurat_object) <- "Group"
Group1_Group2 <- subset(seurat_object, idents = c("Iso Group2", "Group2"))
Group1_Group2$Group <- factor(Group1_Group2$Group, 
                                 levels = c("Group2", "Iso Group2"))
Idents(Group1_Group2) <- "Group"
levels(Group1_Group2)
table(Group1_Group2$bio_replicate, Group1_Group2$Group)
table(Group1_Group2$bio_replicate, Group1_Group2$final_clusters)
deg_group = run_de(Group1_Group2, de_family = 'pseudobulk', de_method = 'edgeR', de_type = 'LRT', 
                   cell_type_col = "final_clusters", replicate_col = "bio_replicate", label_col = "Group")

#quick check on signifcant p values 
num_sig_genes <- deg_group %>%
  filter(p_val_adj < 0.05) %>%
  nrow()

num_sig_genes <- deg_group %>%
  filter(p_val < 0.05) %>%
  nrow()

#save DEGs
write.csv(deg_group,
          "~/directory/deg_group.csv")

#Loop over all clusters for generating DEG output------------------------------
#load run_de output 
deg_group <- read_csv("~/directory/deg_group.csv", col_types = cols(...1 = col_skip()))

parent_directory <- "~/directory"

folder_names <- unique(deg_group$cell_type)
folder_names

# Create subfolders based on the vector
for (folder_name in folder_names) {
  subfolder_path <- file.path(parent_directory, folder_name)
  dir.create(subfolder_path, recursive = TRUE, showWarnings = T)
  if (file.exists(subfolder_path)) {
    cat("Subfolder", folder_name, "created successfully.\n")
  } else {
    warning("Failed to create subfolder", folder_name)
  }
}

# Loop to get upregulated genes in each cluster 
for (cluster_id in unique(deg_group$cell_type)) {
  
  current_folder <- file.path(parent_directory, as.character(cluster_id))

  cluster_df <- deg_group[deg_group$cell_type == cluster_id, ]
  
  folder_name_markers <- as_tibble(cluster_df)

  result_file <- file.path(current_folder, paste0(cluster_id, "_DEGs.csv"))
  write.csv(cluster_df, file = result_file, row.names = FALSE)
 
  folder_name_up_in_Group1 <- subset(folder_name_markers, p_val_adj <= 0.05 & avg_logFC < 0)
  names(folder_name_up_in_Group1)[names(folder_name_up_in_Group1) == "gene"] <- "up_in_Group1"
  write.csv(folder_name_up_in_Group1, row.names = TRUE,
            file.path(current_folder, paste0(cluster_id, "_up_in_Group1.csv")))

  folder_name_up_in_Group2 <- subset(folder_name_markers, p_val_adj <= 0.05 & avg_logFC > 0)
  names(folder_name_up_in_Group2)[names(folder_name_up_in_Group2) == "gene"] <- "up_in_Group2"
  write.csv(folder_name_up_in_Group2, row.names = TRUE,
            file.path(current_folder, paste0(cluster_id, "_up_in_Group2.csv")))
}

#loop to get genes for comparative pw analysis in metascape 
for (folder_name in folder_names) {

  current_folder <- file.path(parent_directory, folder_name)

  up_in_Group1_path <- file.path(current_folder, paste0(folder_name, "_up_in_Group1.csv"))
  up_in_Group2_path <- file.path(current_folder, paste0(folder_name, "_up_in_Group2.csv"))
  
  up_in_Group1 <- read_csv(up_in_Group1_path)
  up_in_Group2 <- read_csv(up_in_Group2_path)
  
  names(up_in_Group1)[names(up_in_Group1) == "gene"] <- "up_in_Group1"
  names(up_in_Group2)[names(up_in_Group2) == "gene"] <- "up_in_Group2"
  
  mscape <- merge(up_in_Group1, up_in_Group2, all = TRUE)

  mscape <- select(mscape, up_in_Group1, up_in_Group2)

  if (!file.exists(current_folder)) {
    dir.create(current_folder, recursive = TRUE)
  }

  excel_path <- file.path(current_folder, paste0(folder_name, "_mscape.xlsx"))
  write.xlsx(mscape, excel_path)
}


#calculate significantly regulated genes for waterfall plot---------------------

#set up function to calculate genes 
calculate_gene_counts <- function(df) {
 
  df <- df[complete.cases(df), ]
  
  upregulated <- df[df$p_val_adj < 0.05 & df$avg_logFC > 0, ]
  downregulated <- df[df$p_val_adj < 0.05 & df$avg_logFC < 0, ]
  not_regulated <- df[df$p_val_adj > 0.05, ]
  
  total_genes <- nrow(df)
  
  num_upregulated <- nrow(upregulated)
  num_downregulated <- nrow(downregulated)
  num_not_regulated <- nrow(not_regulated)

  return(list(
    total_genes = total_genes,
    upregulated = num_upregulated,
    downregulated = num_downregulated,
    not_regulated = num_not_regulated
  ))
}

# Loop over cell types to calculate gene numbers for waterfall plots
for (folder_name in folder_names) {

  subfolder_path <- file.path(parent_directory, folder_name)

  files <- list.files(subfolder_path, pattern = "_DEGs.csv", full.names = TRUE)
  
  for (file in files) {
    
    df <- read.csv(file)

    gene_counts <- calculate_gene_counts(df)

    result_file <- file.path(subfolder_path, paste0("gene_counts_", tools::file_path_sans_ext(basename(file)), ".csv"))
    write.csv(gene_counts, file = result_file)
  }
}

# Read in gene counts for each cluster
gene_counts_list <- list()

# Loop through each cluster to read the gene count CSV file
for (folder_name in folder_names) {
  file_path <- file.path("~/directory", folder_name, paste0("gene_counts_", folder_name, "_DEGs.csv"))
  gene_counts <- read.csv(file_path, row.names = 1) 
  gene_counts_list[[folder_name]] <- gene_counts
}


# Combine gene counts into a single data frame
waterfall_Group1_Group2 <- do.call(rbind, gene_counts_list)
openxlsx::write.xlsx(waterfall_Group1_Group2, asTable =T, rowNames =T,
                     "~/directory/waterfall_Group1_Group2.xlsx")

#check on aggregated counts cluster of interest for vallidation-----------------
library(viridis)
library(ggplot2)
Idents(Group1_Group2) <- "final_clusters"
cluster_1 <- subset(Group1_Group2, idents = "cluster_1")

cluster_1_DEGs <- read_csv("~/directory/cluster_1/cluster_1_DEGs.csv")
cluster_1_DEGs <- subset(cluster_1_DEGs, p_val_adj < 0.05)
cluster_1_DEGs %>%
  dplyr::filter(avg_logFC > 0) %>%
  slice_head(n = 150) %>%
  ungroup() -> top150

cluster_1$Aggregated_ID <- sub("_c|_n", "", cluster_1$bio_replicate)
cluster_1$ID_Group <- paste(cluster_1$bio_replicate, cluster_1$Group, sep = "_")

cluster_averages <- AverageExpression(object = cluster_1,
                                      assay = "RNA3",
                                      return.seurat = TRUE,
                                      group.by = 'ID_Group')

DoHeatmap(cluster_averages,
          features = top150$gene,
          angle = 90,
          raster = FALSE,
          assay = "RNA3",
          slot = "scale.data",
          draw.lines = FALSE,
          label = TRUE,
          lines.width = 1,
          disp.min = -2,
          disp.max = 2,
          size = 5) +
  NoLegend() +
  scale_fill_viridis_c(option = "magma")