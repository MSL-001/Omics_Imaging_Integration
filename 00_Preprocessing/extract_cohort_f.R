args <- commandArgs(trailingOnly = TRUE)

fold_suffix <- args[1]


proteomics <- read.csv("proteomics.csv")
metabolomics <- read.csv("metabolomics.csv")
meanff_80 <- read.csv("female_mean_ff_VIBESegmentator_80.csv")
medianff_80 <- read.csv("female_median_ff_VIBESegmentator_80.csv")
volume_80 <- read.csv("female_vol_cm3_VIBESegmentator_80.csv")
meanff_100 <- read.csv("female_mean_ff_VIBESegmentator_100.csv")
medianff_100 <- read.csv("female_median_ff_VIBESegmentator_100.csv")
volume_100 <- read.csv("female_vol_cm3_VIBESegmentator_100.csv")

meanff_80_processed <- read.csv("female_mean_ff_VIBESegmentator_80_processed.csv")
medianff_80_processed <- read.csv("female_median_ff_VIBESegmentator_80_processed.csv")
volume_80_processed <- read.csv("female_vol_cm3_VIBESegmentator_80_processed.csv")
meanff_100_processed <- read.csv("female_mean_ff_VIBESegmentator_100_processed.csv")
medianff_100_processed <- read.csv("female_median_ff_VIBESegmentator_100_processed.csv")
volume_100_processed <- read.csv("female_vol_cm3_VIBESegmentator_100_processed.csv")

df <- list(meanff_80, medianff_80, volume_80, meanff_100, medianff_100, volume_100, meanff_80_processed,
           medianff_80_processed, volume_80_processed, meanff_100_processed, medianff_100_processed, volume_100_processed)
result <- lapply(df, function(x) {
  colnames(x)[colnames(x) == "subject_ID"] <- "eid"
  x <- x[, names(x) != "unused", drop = FALSE]
  x <- x[, names(x) != "prostate", drop = FALSE]
  return(x)
})

meanff_80   <- result[[1]]
medianff_80 <- result[[2]]
volume_80   <- result[[3]]
meanff_100   <- result[[4]]
medianff_100 <- result[[5]]
volume_100   <- result[[6]]
meanff_80_processed   <- result[[7]]
medianff_80_processed <- result[[8]]
volume_80_processed   <- result[[9]]
meanff_100_processed   <- result[[10]]
medianff_100_processed <- result[[11]]
volume_100_processed   <- result[[12]]

extract_from_df <- function(df, name){
  cohort_df <- merge(df, cohort, by="eid")
  print(paste(name, nrow(cohort_df), ncol(cohort_df), sep=" "))
  write.csv(cohort_df, paste0(name, "_data", fold_suffix,".csv"), row.names=F)
}

cohort <- read.csv("cohort.csv")

extract_from_df(proteomics, "proteomics")
extract_from_df(metabolomics, "metabolomics")
extract_from_df(meanff_80, "meanff_80")
extract_from_df(medianff_80, "medianff_80")
extract_from_df(volume_80, "volume_80")
extract_from_df(meanff_100, "meanff_100")
extract_from_df(medianff_100, "medianff_100")
extract_from_df(volume_100, "volume_100")
extract_from_df(meanff_80_processed , "meanff_80_processed")
extract_from_df(medianff_80_processed , "medianff_80_processed")
extract_from_df(volume_80_processed , "volume_80_processed")
extract_from_df(meanff_100_processed , "meanff_100_processed")
extract_from_df(medianff_100_processed , "medianff_100_processed")
extract_from_df(volume_100_processed , "volume_100_processed")




