library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

train_eids_path <- args[1]
test_eids_path <- args[2]
input_omics_path <- args[3]
input_volume_path <- args[4]
input_ff_path <- args[5]


train_eids  <- read.csv(train_eids_path)
test_eids  <- read.csv(test_eids_path)

omics_data <- read.csv(input_omics_path)
volume_data <- read.csv(input_volume_path)
ff_data <- read.csv(input_ff_path)
age_data <- read.csv("age_data.csv")
height_data <- read.csv("height_data.csv")

colnames(volume_data)[colnames(volume_data) == "subject_ID"] <- "eid"
colnames(ff_data)[colnames(ff_data) == "subject_ID"] <- "eid"
meta_data <- merge(age_data, height_data, by="eid")

meta_train <- meta_data %>%
  filter(eid %in% train_eids$eid)

meta_test <- meta_data %>%
  filter(eid %in% test_eids$eid)

omics_train <- omics_data %>%
  filter(eid %in% train_eids$eid)

omics_test <- omics_data %>%
  filter(eid %in% test_eids$eid)

volume_data_train <- volume_data %>%
  filter(eid %in% train_eids$eid)

volume_data_test <- volume_data %>%
  filter(eid %in% test_eids$eid)

ff_data_train <- ff_data %>%
  filter(eid %in% train_eids$eid)

ff_data_test <- ff_data %>%
  filter(eid %in% test_eids$eid)


meta_train <- meta_train %>% select(-eid)
meta_test <- meta_test %>% select(-eid)
omics_train <- omics_train %>% select(-eid)
omics_test <- omics_test %>% select(-eid)
volume_data_train <- volume_data_train %>% select(-eid)
volume_data_test <- volume_data_test %>% select(-eid)
ff_data_train <- ff_data_train %>% select(-eid)
ff_data_test <- ff_data_test %>% select(-eid)

column_means <- function(df, type){
  cols <- colnames(df)
  result <- data.frame(type=type)
  for (c in cols){
    result[[c]] <- mean(df[[c]], na.rm = TRUE)
  }
  return(result)
}

column_sd <- function(df, type){
  cols <- colnames(df)
  result <- data.frame(type=type)
  for (c in cols){
    result[[c]] <- sd(df[[c]], na.rm = TRUE)
  }
  return(result)
}

train_means_meta <- column_means(meta_train,"meta_train")
test_means_meta <- column_means(meta_test,"meta_test")
train_means_omics <- column_means(omics_train,"omic_train")
test_means_omics <- column_means(omics_test,"omic_test")
train_means_volume <- column_means(volume_data_train,"volume_train")
test_means_volume <- column_means(volume_data_test,"volume_test")
train_means_ff <- column_means(ff_data_train,"ff_train")
test_means_ff <- column_means(ff_data_test,"ff_test")


train_sd_meta <- column_sd(meta_train,"meta_train")
test_sd_meta <- column_sd(meta_test,"meta_test")
train_sd_omics <- column_sd(omics_train,"omic_train")
test_sd_omics <- column_sd(omics_test,"omic_test")
train_sd_volume <- column_sd(volume_data_train,"volume_train")
test_sd_volume <- column_sd(volume_data_test,"volume_test")
train_sd_ff <- column_sd(ff_data_train,"ff_train")
test_sd_ff <- column_sd(ff_data_test,"ff_test")

meta_means <- rbind(train_means_meta, test_means_meta)
meta_sd <- rbind(train_sd_meta, test_sd_meta)

omic_means <- rbind(train_means_omics, test_means_omics)
omic_sd <- rbind(train_sd_omics, test_sd_omics)

image_means <- rbind(train_means_volume, test_means_volume, train_means_ff, test_means_ff)
image_sd <- rbind(train_sd_volume, test_sd_volume, train_sd_ff, test_sd_ff)

write.csv(meta_means, "meta_means.csv", row.names=FALSE)
write.csv(meta_sd, "meta_sd.csv", row.names=FALSE)
write.csv(omic_means, "omic_means.csv", row.names=FALSE)
write.csv(omic_sd, "omic_sd.csv", row.names=FALSE)
write.csv(image_means, "image_means.csv", row.names=FALSE)
write.csv(image_sd, "image_sd.csv", row.names=FALSE)



