library(readxl)
library(pheatmap)
library(dplyr)
library(tidyr)

generate_heatmap <- function(title, cohort, file, generate_cluster){
  data_folder <- "Data/01_Statistics/"

  file_path <- paste0(data_folder, cohort, file)

  file_path_encoding <- paste0("Data/", "Metadata/encoding.csv")

  encoding <- read.csv(file_path_encoding)

  results <- read.csv(file_path)

  #significant_results <- subset(results, Rows_used > 500 &  q_value < 0.01)

  #filtered_results <- subset(significant_results, abs(Adjusted_corr) > 0)

  #omic_keep <- unique(filtered_results$Omics_Feature)
  #feature_keep <- unique(filtered_results$BC_feature)

  #significant_results <- significant_results %>%
    #dplyr::filter(Omics_Feature %in% omic_keep)

  #significant_results <- significant_results %>%
  #  dplyr::filter(BC_feature %in% feature_keep)

  results$Adjusted_corr[!(results$Rows_used > 500 & results$q_value < 0.01)] <- NA

  mat_df <- results %>%
    dplyr::select(Omics_Feature, BC_feature, Adjusted_corr) %>%
    tidyr::pivot_wider(names_from = BC_feature, values_from = Adjusted_corr)

  lookup <- setNames(encoding$Name, encoding$...1)


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
    row_hc <- hclust(dist(mat_cluster))
  }
  #col_hc <- hclust(dist(t(mat_cluster)))

  # plot using original matrix, but force clustering order
  pheatmap(mat,
           cluster_rows = row_hc,
           cluster_cols = F,
           na_col = "grey90",
           fontsize_row= 5,
           main=title)
}

generate_heatmap("Metabolomics Male F1 medianff vibe100", "Met_m/fold1/", "corr_adj_medianff_100.csv", generate_cluster = T)
generate_heatmap("Metabolomics Male F2 medianff vibe100", "Met_m/fold2/", "corr_adj_medianff_100.csv", generate_cluster = F)
generate_heatmap("Metabolomics Male F1 medianff vibe100 processed", "Met_m/fold1/", "corr_adj_medianff_100_processed.csv", generate_cluster = F)
generate_heatmap("Metabolomics Male F2 medianff vibe100 processed", "Met_m/fold2/", "corr_adj_medianff_100_processed.csv", generate_cluster = F)
generate_heatmap("Metabolomics Male F1 meanff vibe100", "Met_m/fold1/", "corr_adj_meanff_100.csv", generate_cluster = F)
generate_heatmap("Metabolomics Male F2 meanff vibe100", "Met_m/fold2/", "corr_adj_meanff_100.csv", generate_cluster = F)
generate_heatmap("Metabolomics Male F1 meanff vibe100", "Met_m/fold1/", "corr_adj_meanff_100_processed.csv", generate_cluster = F)
generate_heatmap("Metabolomics Male F2 meanff vibe100", "Met_m/fold2/", "corr_adj_meanff_100_processed.csv", generate_cluster = F)
generate_heatmap("Metabolomics Male F1 medianff vibe80", "Met_m/fold1/", "corr_adj_medianff_80.csv", generate_cluster = F)
generate_heatmap("Metabolomics Male F2 medianff vibe80", "Met_m/fold2/", "corr_adj_medianff_80.csv", generate_cluster = F)
generate_heatmap("Metabolomics Male F1 medianff vibe80", "Met_m/fold1/", "corr_adj_medianff_80_processed.csv", generate_cluster = F)
generate_heatmap("Metabolomics Male F2 medianff vibe80", "Met_m/fold2/", "corr_adj_medianff_80_processed.csv", generate_cluster = F)
generate_heatmap("Metabolomics Male F1 meanff vibe80", "Met_m/fold1/", "corr_adj_meanff_80.csv", generate_cluster = F)
generate_heatmap("Metabolomics Male F2 meanff vibe80", "Met_m/fold2/", "corr_adj_meanff_80.csv", generate_cluster = F)
generate_heatmap("Metabolomics Male F1 medianff vibe80", "Met_m/fold1/", "corr_adj_medianff_80_processed.csv", generate_cluster = F)
generate_heatmap("Metabolomics Male F2 medianff vibe80", "Met_m/fold2/", "corr_adj_medianff_80_processed.csv", generate_cluster = F)

generate_heatmap("Proteomics Male F1 medianff vibe100", "Single_prot_m/fold1/", "corr_adj_medianff_100.csv", generate_cluster = T)
generate_heatmap("Proteomics Male F2 medianff vibe100", "Single_prot_m/fold2/", "corr_adj_medianff_100.csv", generate_cluster = F)
generate_heatmap("Proteomics Male F1 medianff vibe100 processed", "Single_prot_m/fold1/", "corr_adj_medianff_100_processed.csv", generate_cluster = F)
generate_heatmap("Proteomics Male F2 medianff vibe100 processed", "Single_prot_m/fold2/", "corr_adj_medianff_100_processed.csv", generate_cluster = F)
generate_heatmap("Proteomics Male F1 meanff vibe100", "Single_prot_m/fold1/", "corr_adj_meanff_100.csv", generate_cluster = F)
generate_heatmap("Proteomics Male F2 meanff vibe100", "Single_prot_m/fold2/", "corr_adj_meanff_100.csv", generate_cluster = F)
generate_heatmap("Proteomics Male F1 meanff vibe100", "Single_prot_m/fold1/", "corr_adj_meanff_100_processed.csv", generate_cluster = F)
generate_heatmap("Proteomics Male F2 meanff vibe100", "Single_prot_m/fold2/", "corr_adj_meanff_100_processed.csv", generate_cluster = F)
generate_heatmap("Proteomics Male F1 medianff vibe80", "Single_prot_m/fold1/", "corr_adj_medianff_80.csv", generate_cluster = F)
generate_heatmap("Proteomics Male F2 medianff vibe80", "Single_prot_m/fold2/", "corr_adj_medianff_80.csv", generate_cluster = F)
generate_heatmap("Proteomics Male F1 medianff vibe80", "Single_prot_m/fold1/", "corr_adj_medianff_80_processed.csv", generate_cluster = F)
generate_heatmap("Proteomics Male F2 medianff vibe80", "Single_prot_m/fold2/", "corr_adj_medianff_80_processed.csv", generate_cluster = F)
generate_heatmap("Proteomics Male F1 meanff vibe80", "Single_prot_m/fold1/", "corr_adj_meanff_80.csv", generate_cluster = F)
generate_heatmap("Proteomics Male F2 meanff vibe80", "Single_prot_m/fold2/", "corr_adj_meanff_80.csv", generate_cluster = F)
generate_heatmap("Proteomics Male F1 medianff vibe80", "Single_prot_m/fold1/", "corr_adj_medianff_80_processed.csv", generate_cluster = F)
generate_heatmap("Proteomics Male F2 medianff vibe80", "Single_prot_m/fold2/", "corr_adj_medianff_80_processed.csv", generate_cluster = F)