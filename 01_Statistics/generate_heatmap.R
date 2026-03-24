library(readxl)
library(pheatmap)
library(dplyr)
library(tidyr)

generate_heatmap <- function(title, cohort, file, corr_type, generate_cluster=F, scale=F, rows=F,remove_prostate=F,panel=F){
  data_folder <- "Data/01_Statistics/"

  file_path <- paste0(data_folder, cohort, file)

  file_path_encoding <- paste0("Data/", "Metadata/encoding.csv")

  encoding <- read.csv(file_path_encoding)

  results <- read.csv(file_path)

  Adj_corr <- paste0("Adjusted_corr_", corr_type)
  q_value_col <- paste0("q_value_", corr_type)

  if (!identical(panel, F)) {
    # keep_panels <- c(panel, paste(panel, "II"))
    keep_panels <- c(panel)
    keep_codes <- encoding$Code[encoding$Panel %in% keep_panels]
    results <- results %>% dplyr::filter(Omics_Feature %in% keep_codes)
  }

  # significant_results <- subset(results, Rows_used > 500 &  q_value_col < 0.01)
  #
  # filtered_results <- subset(significant_results, abs(Adj_corr) > 0)
  #
  # omic_keep <- unique(filtered_results$Omics_Feature)
  # feature_keep <- unique(filtered_results$BC_feature)
  #
  # results <- results %>%
  #   dplyr::filter(Omics_Feature %in% omic_keep)

  results[[Adj_corr]][!(results$Rows_used > 1000 & results[[q_value_col]] < 0.01)] <- NA

  mat_df <- results %>%
    dplyr::select(Omics_Feature, BC_feature, all_of(Adj_corr)) %>%
    tidyr::pivot_wider(names_from = BC_feature, values_from = all_of(Adj_corr))

  if (!"cardiac_fat" %in% colnames(mat_df)) {
  mat_df[["cardiac_fat"]] <- NA
  }

  mat_df <- mat_df[, names(mat_df) != "unused", drop = FALSE]

  cols <- colnames(mat_df)
  cols <- append(cols[cols != "cardiac_fat"],
                    "cardiac_fat",
                    after = match("bone_other", cols[cols != "cardiac_fat"]))
  mat_df <- mat_df[, cols]

  if (remove_prostate == T){
    mat_df <- mat_df[, colnames(mat_df) != "prostate"]
  }
  lookup <- setNames(encoding$Name, encoding$Code)


  mat_df$Omics_Feature <- ifelse(
    mat_df$Omics_Feature %in% names(lookup),
    lookup[mat_df$Omics_Feature],
    mat_df$Omics_Feature
  )

  mat <- as.matrix(mat_df[,-1])



  rownames(mat) <- make.unique(mat_df$Omics_Feature)

  mat_cluster <- mat
  mat_cluster[is.na(mat_cluster)] <- 0

  # compute clustering using the NA->0 matrix
  if (generate_cluster==T){
    row_hc <<- hclust(dist(mat_cluster))
  }

  mat <- mat[row_hc$order, , drop = FALSE]
  #col_hc <- hclust(dist(t(mat_cluster)))
  fontsize <- 5
  if (!identical(rows, F)){
    lower <- rows[1]
    higher <- rows[2]
    mat <- mat[lower:higher, , drop = FALSE]
    fontsize <- 10
  }

  colors <- colorRampPalette(c("blue", "white", "red"))(100)
  breaks <- seq(scale[1], scale[2], length.out = 101)

  show_row_names <- nrow(mat) <= 300
  pheatmap(mat,
           cluster_rows = F,
           cluster_cols = F,
           na_col = "grey90",
           breaks = breaks,
           color = colors,
           fontsize_row= fontsize,
           main=title,
           show_rownames = show_row_names)
}

rows <- F
scale <- c(-0.4, 0.4)
cohort <- "Met_m/"

for (feature_type in c("robustmeanff_100_processed", "vol_100_processed")){
  clusters <- T
  for (corr_type in c("pearson_height", "pearson", "spearman", "spearman_height")){
      generate_heatmap( paste0("Metabolomics Male ", feature_type,"_",corr_type), cohort, paste0("corr_adj_",feature_type,".csv"), corr_type, scale=scale,  generate_cluster = clusters, rows=rows)
        clusters <- F
  }
}

cohort <- "Met_f/"
for (feature_type in c("robustmeanff_100_processed", "vol_100_processed")){
  clusters <- T
  for (corr_type in c("pearson_height", "pearson", "spearman", "spearman_height")){
      generate_heatmap( paste0("Metabolomics Female ", feature_type,"_",corr_type), cohort, paste0("corr_adj_",feature_type,".csv"), corr_type, scale=scale,  generate_cluster = clusters, rows=rows)
        clusters <- F
  }
}
scale <- c(-0.4, 0.4)
# for (panel in c("Oncology", "Neurology", "Cardiometabolic", "Inflammation","Oncology II", "Neurology II", "Cardiometabolic II", "Inflammation II")){
#   cohort <- "Single_prot_m/"
#   generate_heatmap(paste("Proteomics Male robust meanff vibe100_processed", panel), cohort, "corr_adj_robustmeanff_100_processed.csv", scale, generate_cluster = T, rows=rows, panel=panel)
#   generate_heatmap(paste("Proteomics Male volume vibe100_processed", panel), cohort, "corr_adj_vol_100_processed.csv", scale, generate_cluster = T, rows=rows, panel=panel)
# }
cohort <- "Single_prot_m/"

for (feature_type in c("robustmeanff_100_processed", "vol_100_processed")){
  clusters <- T
  for (corr_type in c("pearson_height", "pearson", "spearman", "spearman_height")){
      generate_heatmap(paste0("Proteomics Male ", feature_type,"_",corr_type), cohort, paste0("corr_adj_",feature_type,".csv"), corr_type, scale=scale,  generate_cluster = clusters, rows=rows)
        clusters <- F
  }
}


# for (panel in c("Oncology", "Neurology", "Cardiometabolic", "Inflammation","Oncology II", "Neurology II", "Cardiometabolic II", "Inflammation II")){
#   cohort <- "Single_prot_f/"
#   generate_heatmap(paste("Proteomics Female robust meanff vibe100_processed", panel), cohort, "corr_adj_robustmeanff_100_processed.csv", scale, generate_cluster = T, rows=rows, panel=panel)
#   generate_heatmap(paste("Proteomics Female volume vibe100_processed", panel), cohort, "corr_adj_vol_100_processed.csv", scale, generate_cluster = T, rows=rows, panel=panel)
# }

cohort <- "Single_prot_f/"

for (feature_type in c("robustmeanff_100_processed", "vol_100_processed")){
  clusters <- T
  for (corr_type in c("pearson_height", "pearson", "spearman", "spearman_height")){
      generate_heatmap(paste0("Proteomics Female ", feature_type,"_",corr_type), cohort, paste0("corr_adj_",feature_type,".csv"), corr_type, scale=scale,  generate_cluster = clusters, rows=rows)
        clusters <- F
  }
}