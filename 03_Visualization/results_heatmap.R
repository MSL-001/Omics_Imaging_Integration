
library(pheatmap)
library(dplyr)
library(tidyr)
library(readxl)

file_path <- "C:/Users/matti/OneDrive/Documents/Combined Results Heldout.xlsx"
met_m_data <- read_excel(file_path, sheet="Metabolites Male")
met_f_data <- read_excel(file_path, sheet="Metabolites Female")
prot_f_data <- read_excel(file_path, sheet="Proteins Female")
prot_m_data <- read_excel(file_path, sheet="Proteins Male")

met_m_data_use <- data.frame(
  feature = met_m_data[["feature"]],
  met_m_ff = met_m_data[["R2 FF"]],
  met_m_ff_height = met_m_data[["R2 FF Height"]],
  met_m_vol = met_m_data[["R2 Volume"]],
  met_m_vol_height = met_m_data[["R2 Volume Height"]]
)

met_f_data_use <- data.frame(
  feature = met_f_data[["feature"]],
  met_f_ff = met_f_data[["R2 FF"]],
  met_f_ff_height = met_f_data[["R2 FF Height"]],
  met_f_vol = met_f_data[["R2 Volume"]],
  met_f_vol_height = met_f_data[["R2 Volume Height"]]
)

prot_m_data_use <- data.frame(
  feature = prot_m_data[["feature"]],
  prot_m_ff = prot_m_data[["R2 FF"]],
  prot_m_ff_height = prot_m_data[["R2 FF Height"]],
  prot_m_vol = prot_m_data[["R2 Volume"]],
  prot_m_vol_height = prot_m_data[["R2 Volume Height"]]
)

prot_f_data_use <- data.frame(
  feature = prot_f_data[["feature"]],
  prot_f_ff = prot_f_data[["R2 FF"]],
  prot_f_ff_height = prot_f_data[["R2 FF Height"]],
  prot_f_vol = prot_f_data[["R2 Volume"]],
  prot_f_vol_height = prot_f_data[["R2 Volume Height"]]
)

all_results <- merge(met_m_data_use, prot_m_data_use, by="feature")
all_results <- merge(all_results, prot_f_data_use, by="feature", all.x = TRUE)
all_results <- merge(all_results, met_f_data_use, by="feature", all.x = TRUE)
row.names(all_results) <- all_results[["feature"]]

all_results[all_results<(0.0)] <- 0

all_results[all_results>(5.0)] <- 5
all_results["feature"] <- row.names(all_results)

all_results_save <- all_results
all_results <- all_results[, setdiff(names(all_results), "feature")]
all_results <- as.matrix(all_results)

all_results[all_results<(-0.05)] <- -0.05
# all_results_save <- all_results_save[, setdiff(names(all_results_save), setdiff(names(all_results_save),c("feature", "met_m_ff")))]

write.csv(all_results_save, "results for heatmap.csv", row.names=F)
# pheatmap(all_results,
#            cluster_rows = T,
#            cluster_cols = T,
#            na_col = "grey90")