library(readxl)
if(!require("pheatmap")){
  install.packages("pheatmap", repos='http://cran.us.r-project.org')
  library(pheatmap)
}
library(dplyr)
library(tidyr)

protein_data <- read.csv("prot_all.csv")
metabolomic_data <- read.csv("metabolomics.csv")

print("merge")

combined_data <- merge(protein_data, metabolomic_data, by="eid")

NA_heatmap <- function(df, name){
  print(paste0(name, ":start"))
  mat <- as.matrix(df)
  mat[!is.na(mat)] <- 1
  mat_cluster <- mat
  mat_cluster[is.na(mat_cluster)] <- 0

  print(paste0(name, ": kmeans"))
  km <- kmeans(mat_cluster, centers = 6, iter.max = 50, nstart = 10)
  ord_row <- order(km$cluster)
  km2 <- kmeans(t(mat_cluster), centers = 6, iter.max = 50, nstart = 10)
  ord_col <- order(km2$cluster)
  mat <- mat[ord_row, ord_col, drop = FALSE]
  mat_cluster <- mat_cluster[ord_row, ord_col, drop = FALSE]



  print(paste0(name, ":cluster"))
  #row_hc <- hclust(dist(mat_cluster))
  #col_hc <- hclust(dist(t(mat_cluster)))

  print(paste0(name, ":heatmap"))
  pheatmap(mat_cluster,
           color = c("black", "white"),
           breaks = c(-0.5, 0.001, 1.5),
           cluster_rows = F,
           cluster_cols = F,
           show_rownames = F,
           show_colnames = F,
           na_col = "grey90",
           fontsize_row= 5,
           filename=paste0(name, ".png")
  )
}

NA_heatmap(protein_data, "proteomic_NA")
NA_heatmap(metabolomic_data, "metabolomic_NA")
NA_heatmap(combined_data, "all_NA")