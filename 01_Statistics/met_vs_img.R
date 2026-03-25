library(dplyr)
# BiocManager::install("impute", update = FALSE, ask = FALSE)
# library(impute)
BiocManager::install(version = "3.20", update = FALSE, ask = FALSE)
BiocManager::install("mixOmics", update = FALSE, ask = FALSE)
library(mixOmics) # import the mixOmics library

img_type_string <- "test"
interested_feature <- "meanff_1"
p_threshold <- 0.05
correction_formula <- ~ age_0 + age_between

args <- commandArgs(trailingOnly = TRUE)

input_metabolomics_path <- args[1]
input_image_path <- args[2]
train_eids_path <- args[3]
test_eids_path <- args[4]
img_type_string <- args[5]
interested_feature <- args[6]
height_included <- args[7]

if (height_included == TRUE){
  correction_formula <- ~ age_0 + age_between + p50_i0
}

train_eids  <- read.csv(train_eids_path)
test_eids  <- read.csv(test_eids_path)
metabolomics_data <- read.csv(input_metabolomics_path)
image_data <- read.csv(input_image_path)
age_data <- read.csv("age_data.csv")
height_data <- read.csv("height_data.csv")
encoding <- read.csv("encoding.csv")


NA_threshold <- 0.3

colnames(image_data)[colnames(image_data) == "subject_ID"] <- "eid"

lookup <- setNames(encoding$Name, encoding$Code)

colnames(metabolomics_data) <- ifelse(
  colnames(metabolomics_data) %in% names(lookup),
  lookup[colnames(metabolomics_data)],
  colnames(metabolomics_data)
)

metabolomics_train <- metabolomics_data %>%
  filter(eid %in% train_eids$eid)

metabolomics_test <- metabolomics_data %>%
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

run_modeling <- function (met, img){
  met <- met[, colMeans(is.na(met)) <= NA_threshold]
  met <- met[rowMeans(is.na(met)) <= NA_threshold,]

  img <- img[,colMeans(is.na(img)) <= NA_threshold]
  img <- img[rowMeans(is.na(img)) <= NA_threshold,]

  age_data <- na.omit(age_data)
  height_data <- na.omit(height_data)
  meta_data <- merge(age_data, height_data, by="eid")

  keep_eids <- Reduce(intersect, list(met$eid, img$eid, meta_data$eid))

  met <- met[match(keep_eids, met$eid), ]
  img <- img[match(keep_eids, img$eid), ]
  meta_data <- meta_data[match(keep_eids, meta_data$eid), ]


  met_eids <- met$eid

  meta_met_aligned <- meta_data[match(met$eid, meta_data$eid), ]

  met <- met[, setdiff(names(met), "eid"), drop = FALSE]

  met <- as.data.frame(lapply(met, function(v) {
    residualize_vector(v, meta_met_aligned, correction_formula)
  }))

  met$eid <- met_eids

  img_eids <- img$eid

  img <- img[, interested_feature, drop = FALSE]

  meta_img_aligned <- meta_data[match(img_eids, meta_data$eid), ]

  img[[interested_feature]] <- residualize_vector(
    img[[interested_feature]],
    meta_img_aligned,
    correction_formula
  )

  img[["eid"]] <- img_eids

  norm_met <- fit_normalizer(met, "eid")
  met <- apply_normalizer(met, norm_met)

  norm_img <- fit_normalizer(img, "eid")
  img <- apply_normalizer(img, norm_img)

  X <- met[, setdiff(names(met), "eid"), drop = FALSE]
  y <- img[[interested_feature]]

  univar_stats <- data.frame(
    metabolite = colnames(X),
    cor = NA_real_,
    pvalue = NA_real_
  )

  for (j in seq_along(X)) {
    test_X <- X[[j]]
    test_y <- y
    ok <- complete.cases(test_X, test_y)

    test <- suppressWarnings(cor.test(test_X[ok], test_y[ok], method = "pearson"))
    univar_stats$cor[j] <- unname(test$estimate)
    univar_stats$pvalue[j] <- test$p.value
  }

  met[is.na(met)] <- 0
  img <- na.omit(img)

  keep_eids <- Reduce(intersect, list(met$eid, img$eid, meta_data$eid))

  met <- met[match(keep_eids, met$eid), ]
  img <- img[match(keep_eids, img$eid), ]

  X <- met[, setdiff(names(met), "eid"), drop = FALSE]
  row.names(X) <- met$eid

  y <- img[[interested_feature]]

  keep_metabolites <- univar_stats$metabolite[univar_stats$pvalue < p_threshold]

  X <- X[, keep_metabolites, drop = FALSE]

  print(paste("Number of metabolites selected:", ncol(X)))

  pls.result <- pls(X, y, ncomp= 10, mode = 'regression')

  saveRDS(pls.result, paste0(interested_feature, "_results.rds"))

  png(filename=paste0(interested_feature, "_sample_spread.png"))
  plotIndiv(pls.result, ind.names = FALSE)
  dev.off()

  png(filename=paste0(interested_feature, "_Variance.png"))
  plotVar(pls.result, var.names = FALSE)
  dev.off()

  png(filename=paste0(interested_feature, "_Loadings.png"))
  plotLoadings(pls.result, comp = 1, contrib = 'max', method = 'median', block = "X", title="Top 10 Metabolites",ndisplay =10)
  dev.off()

  Q2.pls.result <- perf(pls.result, validation = 'Mfold',
                folds = 10, nrepeat = 5)

  saveRDS(Q2.pls.result, paste0(interested_feature, "_R2_results.rds"))

  png(filename=paste0(interested_feature, "_R2_train.png"))
  vals <- Q2.pls.result$measures$R2$summary$mean
  p <- plot(Q2.pls.result, criterion = "R2", ylim=c(1,-1),title=paste0("R2 of comp 1 is ", round(Q2.pls.result$measures$R2$summary$mean[[1]], 4)))
  p <- p + coord_cartesian(ylim = c(min(min(vals[vals > -1], 0)), max(Q2.pls.result$measures$R2$summary$mean)))
  print(p)
  dev.off()

  png(filename=paste0(interested_feature, "_R2_heldout.png"))
  Q2.pls.result$measures$R2 <- Q2.pls.result$measures$Q2
  vals <- Q2.pls.result$measures$R2$summary$mean
  p <- plot(Q2.pls.result, criterion = "R2", ylim=c(1,-1),title=paste0("R2 of comp 1 is ", round(Q2.pls.result$measures$R2$summary$mean[[1]], 4)))
  p <- p + coord_cartesian(ylim = c(min(min(vals[vals > -1], 0)), max(Q2.pls.result$measures$R2$summary$mean)))
  print(p)
  dev.off()
}

run_test <- function(met_test, met_train, img_test, img_train, n_comp){
  met_train <- met_train[, colMeans(is.na(met_train)) <= NA_threshold]
  met_train <- met_train[rowMeans(is.na(met_train)) <= NA_threshold,]

  img_train <- img_train[,colMeans(is.na(img_train)) <= NA_threshold]
  img_train <- img_train[rowMeans(is.na(img_train)) <= NA_threshold,]

  norm_met <- fit_normalizer(met_train, "eid")
  met_train <- apply_normalizer(met_train, norm_met)
  met_test <- apply_normalizer(met_test, norm_met)

  norm_img <- fit_normalizer(img_train, "eid")
  img_train <- apply_normalizer(img_train, norm_img)
  img_test <- apply_normalizer(img_test, norm_img)


  met_train[is.na(met_train)] <- 0
  met_test[is.na(met_test)] <- 0
  img_train <- na.omit(img_train)
  img_test <- na.omit(img_test)

  met_ids <- met$eid
  img_eids <- img$eid
  keep_eids <- intersect(met_ids, img_eids)

  met <- met[met$eid %in% keep_eids, ]

  img <- img[img$eid %in% keep_eids, ]


}

run_modeling(metabolomics_train, image_data_train)
