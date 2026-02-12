proteomics <- read.csv("prot_all.csv")
metabolomics <- read.csv("metabolomics.csv")
image_features <- read.csv("sex_bmi_meanvol.csv")
colnames(image_features)[colnames(image_features) == "sbj_id"] <- "eid"
image_features <- image_features[, -which(names(image_features) == "unused")]

cohort <- read.csv("cohort.csv")


cohort_prot <- merge(proteomics, cohort, by="eid")
cohort_met <- merge(metabolomics, cohort, by="eid")
cohort_img <- merge(image_features, cohort, by="eid")


print(paste("Prot", nrow(cohort_prot), ncol(cohort_prot), sep=" "))
print(paste("Met", nrow(cohort_met), ncol(cohort_met), sep=" "))
print(paste("Img", nrow(cohort_img), ncol(cohort_img), sep=" "))


write.csv(cohort_prot, "proteomics_data.csv", row.names=F)
write.csv(cohort_met, "metabolomics_data.csv", row.names=F)
write.csv(cohort_img, "volume_data.csv", row.names=F)


