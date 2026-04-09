library(dplyr)
# BiocManager::install("impute", update = FALSE, ask = FALSE)
# library(impute)
BiocManager::install(version = "3.20", update = FALSE, ask = FALSE)
BiocManager::install("mixOmics", update = FALSE, ask = FALSE)
library(mixOmics) # import the mixOmics library

apply_residualizer <- function(v, covar_data, fitted_obj) {
  df <- covar_data
  df$target <- v

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

img_type_string <- "test"
interested_feature <- "meanff_1"
p_threshold <- 0.05
train_eids_path <- "train_cohort.csv"
test_eids_path <- "test_cohort.csv"
correction_formula <- ~ age_0 + age_between
encoding_path <- "Data/Metadata/encoding.csv"
ncomp <- 2

args <- commandArgs(trailingOnly = TRUE)

input_metabolomics_path <- args[1]
input_image_path <- args[2]
train_eids_path <- args[3]
test_eids_path <- args[4]
img_type_string <- args[5]
interested_feature <- args[6]
height_included <- args[7]
ncomp <- as.integer(args[8])

if (height_included == TRUE){
  correction_formula <- ~ age_0 + age_between + p50_i0
}

train_eids  <- read.csv(train_eids_path)
test_eids  <- read.csv(test_eids_path)
metabolomics_data <- read.csv(input_metabolomics_path)
image_data <- read.csv(input_image_path)
age_data <- read.csv("age_data.csv")
encoding_path <- "encoding.csv"
height_data <- read.csv("height_data.csv")

encoding <- read.csv(encoding_path)

met_resid_models <- readRDS(paste0(interested_feature, "_metabolite_residualizers.rds"))
img_train_fit <- readRDS(paste0(interested_feature, "_image_residualizer.rds"))
norm_met <- readRDS(paste0(interested_feature, "_metabolite_normalizers.rds"))
norm_img <- readRDS(paste0(interested_feature, "_image_normalizer.rds"))

model <- readRDS(paste0(interested_feature, "_results.rds"))

colnames(image_data)[colnames(image_data) == "subject_ID"] <- "eid"

met_test <- metabolomics_data %>%
  filter(eid %in% test_eids$eid)

img_test <- image_data %>%
  filter(eid %in% test_eids$eid)

age_data <- na.omit(age_data)
height_data <- na.omit(height_data)
meta_data <- merge(age_data, height_data, by="eid")

keep_eids <- Reduce(intersect, list(met_test$eid, img_test$eid, meta_data$eid))

met_test <- met_test[match(keep_eids, met_test$eid), ]
img_test <- img_test[match(keep_eids, img_test$eid), ]
meta_data <- meta_data[match(keep_eids, meta_data$eid), ]

met_eids <- met_test$eid
meta_met_aligned <- meta_data[match(met_eids, meta_data$eid), ]

available_metabolites <- intersect(colnames(met_test), names(met_resid_models))
missing_metabolites <- setdiff(colnames(met_test), available_metabolites)

if (length(missing_metabolites) > 0) {
  message("Missing residualizer models for: ", paste(missing_metabolites, collapse = ", "))
}

met_test_mat <- met_test[, available_metabolites, drop = FALSE]
met_test_resid <- met_test_mat

for (j in seq_along(met_test_mat)) {
  met_test_resid[[j]] <- apply_residualizer(
    v = met_test_mat[[j]],
    covar_data = meta_met_aligned,
    fitted_obj = met_resid_models[[colnames(met_test_mat)[j]]]
  )
}

met_test_resid$eid <- met_eids

img_test_eids <- img_test$eid

img_test <- img_test[,interested_feature, drop=FALSE]
meta_img_aligned <- meta_data[match(img_test_eids, meta_data$eid), ]

img_test[[interested_feature]] <- apply_residualizer(
  v = img_test[[interested_feature]],
  covar_data = meta_img_aligned,
  fitted_obj = img_train_fit
)

img_test[["eid"]] <- img_test_eids


met <- apply_normalizer(met_test_resid, norm_met)
img <- apply_normalizer(img_test, norm_img)

met[is.na(met)] <- 0
img <- na.omit(img)

keep_eids <- Reduce(intersect, list(met$eid, img$eid, meta_data$eid))

met <- met[match(keep_eids, met$eid), ]
img <- img[match(keep_eids, img$eid), ]

X <- met[, setdiff(names(met), "eid"), drop = FALSE]

keep_metabolites <- model$names$colnames$X

X <- X[, keep_metabolites, drop = FALSE]

row.names(X) <- met$eid

y_true <- img[[interested_feature]]

pred <- predict(model, newdata=X)

y_est <- pred$predict[,,ncomp]

ss_res <- sum((y_true - y_est)^2)
ss_tot <- sum((y_true - mean(y_true))^2)

r2_test <- if (ss_tot == 0) NA_real_ else 1 - ss_res / ss_tot

print(r2_test)

writeLines(
  paste0("R2_test: ", r2_test),
  con = paste0(interested_feature, "_r2_test.txt")
)