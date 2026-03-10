library(dplyr)
# BiocManager::install("impute", update = FALSE, ask = FALSE)
# library(impute)
BiocManager::install(version = "3.20", update = FALSE, ask = FALSE)
BiocManager::install("mixOmics", update = FALSE, ask = FALSE)
library(mixOmics) # import the mixOmics library

img_type_string <- "test"
interested_feature <- "meanff_1"
p_threshold <- 0.05
train_eids_path <- "train_cohort.csv"
test_eids_path <- "test_cohort.csv"
correction_formula <- ~ age_0 + age_between


args <- commandArgs(trailingOnly = TRUE)

input_proteomics_path <- args[1]
input_image_path <- args[2]
train_eids_path <- args[3]
test_eids_path <- args[4]
img_type_string <- args[5]
interested_feature <- args[6]

train_eids  <- read.csv(train_eids_path)
test_eids  <- read.csv(test_eids_path)
proteomics_data <- read.csv(input_proteomics_path)
image_data <- read.csv(input_image_path)
age_data <- read.csv("age_data.csv")
# height_data <- read.csv("height_data.csv")

NA_threshold <- 0.3

colnames(image_data)[colnames(image_data) == "subject_ID"] <- "eid"

proteomics_train <- proteomics_data %>%
  filter(eid %in% train_eids$eid)

proteomics_test <- proteomics_data %>%
  filter(eid %in% test_eids$eid)

image_data_train <- image_data %>%
  filter(eid %in% train_eids$eid)

image_data_test <- image_data %>%
  filter(eid %in% test_eids$eid)

fit_normalizer <- function(df, keep_out){
  cols <- setdiff(names(df), keep_out)

  mu <- sapply(df[cols], function(x) mean(x, na.rm = TRUE))
  sdv <- sapply(df[cols], function(x) sd(x, na.rm = TRUE))

  sdv[is.na(sdv) | sdv == 0] <- NA_real_

  list(cols = cols, mu = mu, sd = sdv)
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

residualize_vector <- function(v, covar_data, correction_formula) {
  df <- covar_data
  df$target <- v

  fit <- lm(update(correction_formula, target ~ .),
            data = df,
            na.action = na.exclude)

  residuals(fit)
}

run_modeling <- function (prot, img){
  prot <- prot[, colMeans(is.na(prot)) <= NA_threshold]
  prot <- prot[rowMeans(is.na(prot)) <= NA_threshold,]

  img <- img[,colMeans(is.na(img)) <= NA_threshold]
  img <- img[rowMeans(is.na(img)) <= NA_threshold,]

  prot_eids <- prot$eid
  img_eids <- img$eid
  age_data <- na.omit(age_data)
  age_eids <- age_data$eid
  keep_eids <- intersect(prot_eids, img_eids)
  keep_eids <- intersect(keep_eids, age_eids)


  prot <- prot[prot$eid %in% keep_eids, ]

  img <- img[img$eid %in% keep_eids, ]

  age_data <- age_data[age_data$eid %in% keep_eids, ]


  prot_eids <- prot$eid

  age_prot_aligned <- age_data[match(prot$eid, age_data$eid), ]

  prot <- prot[, setdiff(names(prot), "eid"), drop = FALSE]

  prot <- as.data.frame(lapply(prot, function(v) {
    residualize_vector(v, age_prot_aligned, correction_formula)
  }))

  prot$eid <- prot_eids

  img_eids <- img$eid

  img <- img[, interested_feature, drop = FALSE]

  age_img_aligned <- age_data[match(img_eids, age_data$eid), ]

  img[[interested_feature]] <- residualize_vector(
    img[[interested_feature]],
    age_img_aligned,
    correction_formula
  )

  img[["eid"]] <- img_eids

  norm_prot <- fit_normalizer(prot, "eid")
  prot <- apply_normalizer(prot, norm_prot)

  norm_img <- fit_normalizer(img, "eid")
  img <- apply_normalizer(img, norm_img)

  X <- prot[, setdiff(names(prot), "eid"), drop = FALSE]
  y <- img[[interested_feature]]

  univar_stats <- data.frame(
    protein = colnames(X),
    cor = NA_real_,
    pvalue = NA_real_
  )

  for (j in seq_along(X)) {
    test <- suppressWarnings(cor.test(X[[j]], y, method = "pearson"))
    univar_stats$cor[j] <- unname(test$estimate)
    univar_stats$pvalue[j] <- test$p.value
  }

  prot[is.na(prot)] <- 0
  img <- na.omit(img)
  age_data <- na.omit(age_data)

  prot_eids <- prot$eid
  img_eids <- img$eid
  age_eids <- age_data$eid
  keep_eids <- intersect(prot_eids, img_eids)
  keep_eids <- intersect(keep_eids, age_eids)

  prot <- prot[prot$eid %in% keep_eids, ]

  img <- img[img$eid %in% keep_eids, ]

  X <- prot[, setdiff(names(prot), "eid"), drop = FALSE]
  row.names(X) <- prot$eid

  y <- img[[interested_feature]]



  keep_proteins <- univar_stats$protein[univar_stats$pvalue < p_threshold]

  X <- X[, keep_proteins, drop = FALSE]

  print(paste("Number of proteins selected:", ncol(X)))

  pls.result <- pls(X, y, ncomp= 10, mode = 'regression')

  png(filename=paste0("Ind_", img_type_string, ".png"))
  plotIndiv(pls.result, ind.names = FALSE)
  dev.off()

  png(filename=paste0("Var_", img_type_string, ".png"))
  plotVar(pls.result, var.names = FALSE)
  dev.off()

  Q2.pls.result <- perf(pls.result, validation = 'Mfold',
                folds = 10, nrepeat = 5)

  png(filename=paste0("Q2_", img_type_string, ".png"))
  p <- plot(Q2.pls.result, criterion = "Q2")
  print(p)
  dev.off()

  png(filename=paste0("R2_", img_type_string, ".png"))
  p <- plot(Q2.pls.result, criterion = "R2")
  print(p)
  dev.off()
}

run_test <- function(prot_test, prot_train, img_test, img_train, n_comp){
  prot_train <- prot_train[, colMeans(is.na(prot_train)) <= NA_threshold]
  prot_train <- prot_train[rowMeans(is.na(prot_train)) <= NA_threshold,]

  img_train <- img_train[,colMeans(is.na(img_train)) <= NA_threshold]
  img_train <- img_train[rowMeans(is.na(img_train)) <= NA_threshold,]

  norm_prot <- fit_normalizer(prot_train, "eid")
  prot_train <- apply_normalizer(prot_train, norm_prot)
  prot_test <- apply_normalizer(prot_test, norm_prot)

  norm_img <- fit_normalizer(img_train, "eid")
  img_train <- apply_normalizer(img_train, norm_img)
  img_test <- apply_normalizer(img_test, norm_img)


  prot_train[is.na(prot_train)] <- 0
  prot_test[is.na(prot_test)] <- 0
  img_train <- na.omit(img_train)
  img_test <- na.omit(img_test)

  prot_eids <- prot$eid
  img_eids <- img$eid
  keep_eids <- intersect(prot_eids, img_eids)

  prot <- prot[prot$eid %in% keep_eids, ]

  img <- img[img$eid %in% keep_eids, ]


}

run_modeling(proteomics_train, image_data_train)
