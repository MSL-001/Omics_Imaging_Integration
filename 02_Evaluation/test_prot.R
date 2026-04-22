library(dplyr)
# BiocManager::install("impute", update = FALSE, ask = FALSE)
# library(impute)
BiocManager::install(version = "3.20", update = FALSE, ask = FALSE)
BiocManager::install("mixOmics", update = FALSE, ask = FALSE)
library(mixOmics) # import the mixOmics library

apply_residualizer <- function(v, covar_data, fitted_obj) {
  df <- covar_data

  pred <- predict(fitted_obj$fit, newdata = df)
  v - pred
}

apply_normalizer <- function(df, norm){
  cols <- norm$cols
  cols <- intersect(cols, names(df))

  for (c in cols){
    s <- norm$sd[[c]]
    if (is.na(s)) next
    df[[c]] <- (df[[c]] - norm$mu[[c]]) / s
  }
  return(df)
}

args <- commandArgs(trailingOnly = TRUE)

input_proteomics_path <- args[1]
input_image_path <- args[2]
train_eids_path <- args[3]
test_eids_path <- args[4]
feature_type <- args[5]
cohort <- args[6]
gender <- args[7]

features <- c("adrenal_gland_left",
"adrenal_gland_right",
"cardiac_fat",
"duodenum",
"esophagus",
"gallbladder",
"heart",
"inner_fat",
"intestine",
"kidney_left",
"kidney_right",
"liver",
"pancreas",
"sacrum",
"spleen",
"stomach",
"subcutaneous_fat",
"thyroid_gland",
"trachea",
"urinary_bladder",
"lung_upper_lobe_left",
"lung_lower_lobe_left",
"lung_upper_lobe_right",
"lung_middle_lobe_right",
"lung_lower_lobe_right",
"prostate",
"aorta",
"pulmonary_vein",
"brachiocephalic_trunk",
"subclavian_artery_right",
"subclavian_artery_left",
"common_carotid_artery_right",
"common_carotid_artery_left",
"brachiocephalic_vein_left",
"brachiocephalic_vein_right",
"atrial_appendage_left",
"superior_vena_cava",
"inferior_vena_cava",
"portal_vein_and_splenic_vein",
"iliac_artery_left",
"iliac_artery_right",
"iliac_vena_left",
"iliac_vena_right",
"humerus_left",
"humerus_right",
"scapula_left",
"scapula_right",
"clavicula_left",
"clavicula_right",
"femur_left",
"femur_right",
"hip_left",
"hip_right",
"spinal_cord",
"gluteus_maximus_left",
"gluteus_maximus_right",
"gluteus_medius_left",
"gluteus_medius_right",
"gluteus_minimus_left",
"gluteus_minimus_right",
"autochthon_left",
"autochthon_right",
"iliopsoas_left",
"iliopsoas_right",
"sternum",
"costal_cartilages",
"muscle",
"IVD",
"vertebra_body",
"vertebra_posterior_elements",
"spinal_channel",
"bone_other",
"gluteus_maximus_left_fat",
"gluteus_maximus_right_fat",
"gluteus_medius_left_fat",
"gluteus_medius_right_fat",
"gluteus_minimus_left_fat",
"gluteus_minimus_right_fat",
"autochthon_left_fat",
"autochthon_right_fat",
"iliopsoas_left_fat",
"iliopsoas_right_fat",
"muscle_fat")

if (gender == "female"){
  features <- setdiff(features, "prostate")
}

train_eids  <- read.csv(train_eids_path)
test_eids  <- read.csv(test_eids_path)
proteomics_data <- read.csv(input_proteomics_path)
image_data <- read.csv(input_image_path)
age_data <- read.csv("age_data.csv")
encoding_path <- "encoding.csv"
height_data <- read.csv("height_data.csv")

encoding <- read.csv(encoding_path)

colnames(image_data)[colnames(image_data) == "subject_ID"] <- "eid"

omic_test <- proteomics_data %>%
  filter(eid %in% test_eids$eid)

img_test <- image_data %>%
  filter(eid %in% test_eids$eid)

age_data <- na.omit(age_data)
height_data <- na.omit(height_data)
meta_data <- merge(age_data, height_data, by="eid")

meta_data <- meta_data %>%
  filter(eid %in% test_eids$eid)

keep_eids <- Reduce(intersect, list(omic_test$eid, img_test$eid, meta_data$eid))

omic_test_global <- omic_test[match(keep_eids, omic_test$eid), ]
img_test_global <- img_test[match(keep_eids, img_test$eid), ]
meta_data_global <- meta_data[match(keep_eids, meta_data$eid), ]


results <- data.frame(
   feature = features,
   r2_comp1 = NA,
   r2_comp2 = NA,
   r2_comp3 = NA,
   r2_comp4 = NA,
   r2_comp5 = NA,
   r2_comp6 = NA,
   r2_comp7 = NA,
   r2_comp8 = NA,
   r2_comp9 = NA,
   r2_comp10 = NA,
   RMSE_comp1 =NA,
   RMSE_comp2 =NA,
   RMSE_comp3 =NA,
   RMSE_comp4 =NA,
   RMSE_comp5 =NA,
   RMSE_comp6 =NA,
   RMSE_comp7 =NA,
   RMSE_comp8 =NA,
   RMSE_comp9 =NA,
   RMSE_comp10 =NA
)

for (feature in features){
  omic_test <- omic_test_global
  img_test <- img_test_global
  meta_data <- meta_data_global

  omic_resid_models <- readRDS(paste0("/mnt/project/Data/", cohort,"/",feature_type,"/",feature,"/",feature, "_protein_residualizers.rds"))
  img_train_fit <- readRDS(paste0("/mnt/project/Data/", cohort,"/",feature_type,"/",feature,"/",feature, "_image_residualizer.rds"))
  norm_omic <- readRDS(paste0("/mnt/project/Data/", cohort,"/",feature_type,"/",feature,"/",feature, "_protein_normalizers.rds"))
  norm_img <- readRDS(paste0("/mnt/project/Data/", cohort,"/",feature_type,"/",feature,"/",feature, "_image_normalizer.rds"))

  model <- readRDS(paste0("/mnt/project/Data/", cohort,"/",feature_type,"/",feature,"/",feature, "_results.rds"))

  omic_eids <- omic_test$eid
  meta_omic_aligned <- meta_data[match(omic_eids, meta_data$eid), ]

  available_omics <- intersect(colnames(omic_test), names(omic_resid_models))
  missing_omics <- setdiff(colnames(omic_test), available_omics)

  if (length(missing_omics) > 0) {
    message("Missing residualizer models for: ", paste(missing_omics, collapse = ", "))
  }

  omic_test_mat <- omic_test[, available_omics, drop = FALSE]
  omic_test_resid <- omic_test_mat

  for (j in seq_along(omic_test_mat)) {
    omic_test_resid[[j]] <- apply_residualizer(
      v = omic_test_mat[[j]],
      covar_data = meta_omic_aligned,
      fitted_obj = omic_resid_models[[colnames(omic_test_mat)[j]]]
    )
  }

  omic_test_resid$eid <- omic_eids

  img_test_eids <- img_test$eid

  img_test <- img_test[,feature, drop=FALSE]
  meta_img_aligned <- meta_data[match(img_test_eids, meta_data$eid), ]

  img_test[[feature]] <- apply_residualizer(
    v = img_test[[feature]],
    covar_data = meta_img_aligned,
    fitted_obj = img_train_fit
  )

  img_test[["eid"]] <- img_test_eids

  omic <- apply_normalizer(omic_test_resid, norm_omic)
  img <- apply_normalizer(img_test, norm_img)

  omic[is.na(omic)] <- 0
  img <- na.omit(img)

  keep_eids <- Reduce(intersect, list(omic$eid, img$eid, meta_data$eid))

  omic <- omic[match(keep_eids, omic$eid), ]
  img <- img[match(keep_eids, img$eid), ]

  X <- omic[, setdiff(names(omic), "eid"), drop = FALSE]

  keep_omics <- model$names$colnames$X

  X <- X[, keep_omics, drop = FALSE]

  row.names(X) <- omic$eid

  y_true <- img[[feature]]

  pred <- predict(model, newdata=X)

  row_idx <- results$feature == feature

  for (ncomp in c(1,2,3,4,5,6,7,8,9,10)){
    y_est <- pred$predict[,,ncomp]

    rmse <- sqrt(mean((y_true - y_est)^2))

    ss_res <- sum((y_true - y_est)^2)
    ss_tot <- sum(y_true^2)

    r2_test <- if (ss_tot == 0) NA_real_ else 1 - ss_res / ss_tot

    results[row_idx, paste0("r2_comp", ncomp)] <- r2_test
    results[row_idx, paste0("RMSE_comp", ncomp)] <- rmse
  }

}

write.csv(results, paste0(cohort, "_", feature_type, "_test_results.csv"), row.names=FALSE)