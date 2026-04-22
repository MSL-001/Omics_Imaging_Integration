
args <- commandArgs(trailingOnly = TRUE)
feature_type <- args[1]
cohort <- args[2]
gender <- args[3]

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
  feature_training_result <- readRDS(paste0("/mnt/project/Data/", cohort,"/",feature_type,"/",feature,"/",feature,"_R2_results_2.rds"))
  row_idx <- results$feature == feature

  results[row_idx, paste0("r2_comp", 1:10)] <- feature_training_result$summary$q2_mean[1:10]
  results[row_idx, paste0("RMSE_comp", 1:10)] <- feature_training_result$summary$rmse_mean[1:10]
}


write.csv(results, paste0(cohort, "_", feature_type, "_results.csv"), row.names=FALSE)