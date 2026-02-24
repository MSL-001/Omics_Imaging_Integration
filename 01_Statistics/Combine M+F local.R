M_corr <- read.csv("C:/Users/matti/Downloads/correlations_pilot_met_m.csv")

F_corr <- read.csv("C:/Users/matti/Downloads/correlations_pilot_met_f.csv")

library(dplyr)

combined <- M_corr %>%
  left_join(F_corr, by = c("Omics_Feature", "BC_feature", "Feature_type"))

write.csv(combined, "C:/Users/matti/Downloads/correlations_pilot_met_combined.csv", row.names = F)