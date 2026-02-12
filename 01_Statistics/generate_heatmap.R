library(readxl)
library(pheatmap)
library(dplyr)
library(tidyr)

data_folder <- "Data/01_Statistics/"
file <- "m_met_volume.xlsx"
file_path <- paste0(data_folder, file)

encoding <- read_excel(file_path, sheet = "Encoding")


results <- read_excel(file_path, sheet = "Raw_Correlations")

significant_results <- subset(results, Rows_used > 500 &  q_value < 0.01)

filtered_results <- subset(significant_results, abs(Adjusted_corr) > 0)

omic_keep <- unique(filtered_results$Omics_Feature)
feature_keep <- unique(filtered_results$BC_feature)

significant_results <- significant_results %>%
  dplyr::filter(Omics_Feature %in% omic_keep)

significant_results <- significant_results %>%
  dplyr::filter(BC_feature %in% feature_keep)


mat_df <- significant_results %>%
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
row_hc <- hclust(dist(mat_cluster))
col_hc <- hclust(dist(t(mat_cluster)))

# plot using original matrix, but force clustering order
pheatmap(mat,
         cluster_rows = row_hc,
         cluster_cols = col_hc,
         na_col = "grey90",
         fontsize_row= 5)

mat[abs(mat) < 0.2] <- NA
'
pheatmap(mat,
         cluster_rows = row_hc,
         cluster_cols = col_hc,
         na_col = "grey90",
         fontsize_row= 5)
'