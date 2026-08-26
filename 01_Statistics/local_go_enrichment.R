# BiocManager::install("GOfuncR")
# if (!requireNamespace('BiocManager', quietly = TRUE))
# 	install.packages('BiocManager')
# BiocManager::install('Homo.sapiens')

library(GOfuncR)
library(dplyr)

Data_directory <- "Data/01_Statistics/"
cohort <- "Single_prot_f"
feature_type <- "robustmeanff_processed"
# feature_type <- "volume_processed"
feature <- "cardiac_fat"

file_path <- paste0(Data_directory, cohort, "/", feature_type, "/", feature, "_results.rds")

df <- readRDS(file_path)

loadings <- data.frame(
  feature = rownames(df$loadings[["X"]]),
  comp1 = df$loadings[["X"]][, "comp1"]
)

loadings_sorted <- loadings[order(-abs(loadings$comp1)), ]
loadings_positive_rank <- loadings[order(-loadings$comp1), ]
loadings_negative_rank <- loadings[order(loadings$comp1), ]

all_prots <- data.frame(
  Gene_Ids = loadings$feature,
  Is_candidate = 0
)

n_prot <- length(loadings_sorted$feature)

perc_5 <- ceiling(n_prot*0.05)

perc_10 <- ceiling(n_prot*0.1)

perc_15 <- ceiling(n_prot*0.15)

enrichment <- function(sorted_loadings, prots, amount){
  top_prots <- sorted_loadings$feature[1:amount]
  top_prot_input <- prots
  top_prot_input$Is_candidate[top_prot_input$Gene_Ids %in% top_prots] <- 1
  top_prot_input$Gene_Ids <- toupper(top_prot_input$Gene_Ids)
  Go_Enrich_Out<- go_enrich(top_prot_input, silent=T, domains=c("biological_process", "molecular_function"))
  Results<-Go_Enrich_Out$results
  Over_Representation<-Results[Results$FWER_overrep<=0.05,]
  return(Over_Representation)
}

print("Positive 5 percent")
Over_rep_pos_5 <- enrichment(loadings_positive_rank, all_prots, perc_5)
print("Positive 10 percent")
Over_rep_pos_10 <- enrichment(loadings_positive_rank, all_prots, perc_10)
print("Positive 15 percent")
Over_rep_pos_15 <- enrichment(loadings_positive_rank, all_prots, perc_15)

print("Negative 5 percent")
Over_rep_neg_5 <- enrichment(loadings_negative_rank, all_prots, perc_5)
print("Negative 10 percent")
Over_rep_neg_10 <- enrichment(loadings_negative_rank, all_prots, perc_10)
print("Negative 15 percent")
Over_rep_neg_15 <- enrichment(loadings_negative_rank, all_prots, perc_15)

print("All 5 percent")
Over_rep_5 <- enrichment(loadings_sorted, all_prots, perc_5)
print("All 10 percent")
Over_rep_10 <- enrichment(loadings_sorted, all_prots, perc_10)
print("All 15 percent")
Over_rep_15 <- enrichment(loadings_sorted, all_prots, perc_15)

