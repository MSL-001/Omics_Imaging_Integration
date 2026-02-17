library(dplyr)
library(tidyr)


olink_LOD_data <- read.csv("Data/00_Preprocessing/olink_LOD.dat", sep="\t")

olink_LOD_data$Assay <- tolower(olink_LOD_data$Assay)



olink_LOD_data <- olink_LOD_data %>%
  dplyr::filter(Instance == 0)

olink_LOD_data <- olink_LOD_data %>%
  dplyr::select(PlateID, Assay, LOD) %>%
  tidyr::pivot_wider(names_from = Assay, values_from = LOD)


write.csv(olink_LOD_data, "olink_LOD_i0.csv", row.names = F)

olink_LOD_data <- read.csv("Data/00_Preprocessing/olink_LOD.dat", sep="\t")

olink_LOD_data$Assay <- tolower(olink_LOD_data$Assay)

olink_LOD_data <- olink_LOD_data %>%
  dplyr::filter(Instance == 2)

olink_LOD_data <- olink_LOD_data %>%
  dplyr::select(PlateID, Assay, LOD) %>%
  tidyr::pivot_wider(names_from = Assay, values_from = LOD)


write.csv(olink_LOD_data, "olink_LOD_i2.csv", row.names = F)