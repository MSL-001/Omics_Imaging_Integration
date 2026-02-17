proteomics <- read.csv("proteomics.csv")
metabolomics <- read.csv("metabolomics.csv")
volume_80 <- read.csv("sex_bmi_meanvol.csv")
colnames(volume_80)[colnames(volume_80) == "sbj_id"] <- "eid"
meanff_100 <- read.csv("female_mean_ff_VIBESegmentator_100.csv")
medianff_100 <- read.csv("female_median_ff_VIBESegmentator_100.csv")
volume_100 <- read.csv("female_vol_cm3_VIBESegmentator_100.csv")

df <- list(meanff_100, medianff_100, volume_100, volume_80)
result <- lapply(df, function(x) {
  colnames(x)[colnames(x) == "subject_ID"] <- "eid"
  x <- x[, names(x) != "unused", drop = FALSE]
  x <- x[, names(x) != "prostate", drop = FALSE]
  return(x)
})

meanff_100   <- result[[1]]
medianff_100 <- result[[2]]
volume_100   <- result[[3]]
volume_80    <- result[[4]]


volume_80 <- volume_80[, names(volume_80) != "sex", drop = FALSE]
volume_80 <- volume_80[, names(volume_80) != "bmi", drop = FALSE]


cohort <- read.csv("cohort.csv")


cohort_prot <- merge(proteomics, cohort, by="eid")
cohort_met <- merge(metabolomics, cohort, by="eid")
cohort_volume_80 <- merge(volume_80, cohort, by="eid")
cohort_meanff_100 <- merge(meanff_100, cohort, by="eid")
cohort_medianff_100 <- merge(medianff_100, cohort, by="eid")
cohort_volume_100 <- merge(volume_100, cohort, by="eid")



print(paste("Prot", nrow(cohort_prot), ncol(cohort_prot), sep=" "))
print(paste("Met", nrow(cohort_met), ncol(cohort_met), sep=" "))
print(paste("Vol_80", nrow(cohort_volume_80), ncol(cohort_volume_80), sep=" "))
print(paste("Vol_100", nrow(cohort_volume_100), ncol(cohort_volume_100), sep=" "))
print(paste("meanff_100", nrow(cohort_meanff_100), ncol(cohort_meanff_100), sep=" "))
print(paste("medianff_100", nrow(cohort_medianff_100), ncol(cohort_medianff_100), sep=" "))



write.csv(cohort_prot, "proteomics_data.csv", row.names=F)
write.csv(cohort_met, "metabolomics_data.csv", row.names=F)
write.csv(cohort_volume_80, "volume_80_data.csv", row.names=F)
write.csv(cohort_volume_100, "volume_100_data.csv", row.names=F)
write.csv(cohort_meanff_100, "meanff_100_data.csv", row.names=F)
write.csv(cohort_medianff_100, "medianff_100_data.csv", row.names=F)



