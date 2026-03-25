BiocManager::install("GOfuncR")
if (!requireNamespace('BiocManager', quietly = TRUE))
	install.packages('BiocManager')
BiocManager::install('Homo.sapiens')

library(GOfuncR)
library(dplyr)

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

all_prots <- data.frame(
  Gene_Ids = loadings$feature,
  Is_candidate = 0
)

top_150_prots <- loadings_positive_rank$feature[1:150]
top_150_prot_input <- all_prots
top_150_prot_input$Is_candidate[top_150_prot_input$Gene_Ids %in% top_150_prots] <- 1
top_150_prot_input$Gene_Ids <- toupper(top_150_prot_input$Gene_Ids)
Go_Enrich_Out<- go_enrich(top_150_prot_input, domains="biological_process", silent=T)
Results<-Go_Enrich_Out$results
Over_Representation<-Results[Results$FWER_overrep<=0.1,]

top_200_prots <- loadings_positive_rank$feature[1:200]
top_200_prot_input <- all_prots
top_200_prot_input$Is_candidate[top_200_prot_input$Gene_Ids %in% top_200_prots] <- 1
top_200_prot_input$Gene_Ids <- toupper(top_200_prot_input$Gene_Ids)
Go_Enrich_Out_200<- go_enrich(top_200_prot_input, domains="biological_process", silent=T)
Results_200<-Go_Enrich_Out_200$results
Over_Representation_200<-Results_200[Results_200$FWER_overrep<=0.1,]

top_150_prots <- loadings_negative_rank$feature[1:150]
top_150_prot_input <- all_prots
top_150_prot_input$Is_candidate[top_150_prot_input$Gene_Ids %in% top_150_prots] <- 1
top_150_prot_input$Gene_Ids <- toupper(top_150_prot_input$Gene_Ids)
Go_Enrich_Out_neg<- go_enrich(top_150_prot_input, domains="biological_process", silent=T)
Results_neg<-Go_Enrich_Out_neg$results
Over_Representation_neg<-Results_neg[Results_neg$FWER_overrep<=0.1,]

top_200_prots <- loadings_negative_rank$feature[1:200]
top_200_prot_input <- all_prots
top_200_prot_input$Is_candidate[top_200_prot_input$Gene_Ids %in% top_200_prots] <- 1
top_200_prot_input$Gene_Ids <- toupper(top_200_prot_input$Gene_Ids)
Go_Enrich_Out_200_neg<- go_enrich(top_200_prot_input, domains="biological_process", silent=T)
Results_200_neg<-Go_Enrich_Out_200_neg$results
Over_Representation_200_neg<-Results_200_neg[Results_200_neg$FWER_overrep<=0.1,]
