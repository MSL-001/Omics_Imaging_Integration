eids_plates <- read.csv("proteomic_platesi2.csv")

colnames(eids_plates)[colnames(eids_plates) == "p30901_i2"] <- "PlateID"

eids_plates <- na.omit(eids_plates)

plates_LOD <- read.csv("olink_LOD_i2.csv")

merged <- merge(eids_plates, plates_LOD, all.x=T)

print(dim(eids_plates))
print(dim(merged))

merged <- merged[, names(merged) != "PlateID", drop = FALSE]

write.csv(merged, "eids_LOD_i2.csv", row.names = F)

