Data_directory <- "Data/01_Statistics/"
cohort <- "Single_prot_f"
feature_type <- "volume_processed"
feature <- "liver"

file_path <- paste0(Data_directory, cohort, "/", feature_type, "/", feature, "_results.rds")
out_file <- paste0(Data_directory, cohort, "/", feature_type, "/", feature, "_top_loadings.csv")


df <- readRDS(file_path)

loadings <- data.frame(
  feature = rownames(df$loadings[["X"]]),
  comp1 = df$loadings[["X"]][, "comp1"]
)

loadings_sorted <- loadings[order(-abs(loadings$comp1)), ]
loadings_positive_rank <- loadings[order(-loadings$comp1), ]
loadings_negative_rank <- loadings[order(loadings$comp1), ]

write.csv(loadings_sorted, out_file, row.names = F)

