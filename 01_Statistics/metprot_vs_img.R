library(dplyr)
# BiocManager::install("impute", update = FALSE, ask = FALSE)
# library(impute)
BiocManager::install(version = "3.20", update = FALSE, ask = FALSE)
BiocManager::install("mixOmics", update = FALSE, ask = FALSE)
library(mixOmics)

img_type_string <- "test"
interested_feature <- "meanff_1"
p_threshold <- 0.05
train_eids_path <- "train_cohort.csv"
test_eids_path <- "test_cohort.csv"
correction_formula <- ~ age_0 + age_between
encoding_path <- "Data/Metadata/encoding.csv"

args <- commandArgs(trailingOnly = TRUE)

input_proteomics_path <- args[1]
input_metabolomics_path <- args[2]
input_image_path <- args[3]
train_eids_path <- args[4]
test_eids_path <- args[5]
img_type_string <- args[6]
interested_feature <- args[7]
height_included <- args[8]

if (height_included == TRUE){
  correction_formula <- ~ age_0 + age_between + p50_i0
}

train_eids  <- read.csv(train_eids_path)
test_eids  <- read.csv(test_eids_path)
proteomics_data <- read.csv(input_proteomics_path)
metabolomics_data <- read.csv(input_metabolomics_path)
image_data <- read.csv(input_image_path)
age_data <- read.csv("age_data.csv")
encoding_path <- "encoding.csv"
height_data <- read.csv("height_data.csv")

encoding <- read.csv(encoding_path)




NA_threshold <- 0.3

lookup <- setNames(encoding$Name, encoding$Code)

# colnames(proteomics_data) <- ifelse(
#   colnames(proteomics_data) %in% names(lookup),
#   lookup[colnames(proteomics_data)],
#   colnames(proteomics_data)
# )

colnames(metabolomics_data) <- ifelse(
  colnames(metabolomics_data) %in% names(lookup),
  lookup[colnames(metabolomics_data)],
  colnames(metabolomics_data)
)

colnames(image_data)[colnames(image_data) == "subject_ID"] <- "eid"

metabolomics_train <- metabolomics_data %>%
  filter(eid %in% train_eids$eid)

metabolomics_test <- metabolomics_data %>%
  filter(eid %in% test_eids$eid)

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

cv_block_pls_ncomp <- function(X, y, folds = 10, nrepeat = 5, ncomp_max = NULL, seed = 123) {
  ref_ids <- rownames(y)
  if (is.null(ref_ids)) stop("y must have rownames")

  for (b in names(X)) {
    if (is.null(rownames(X[[b]]))) {
      stop(paste("Block", b, "has no rownames"))
    }
    if (!identical(rownames(X[[b]]), ref_ids)) {
      stop(paste("Rownames mismatch between y and block", b))
    }
  }

  set.seed(seed)

  all_results <- list()
  res_i <- 1

  n <- nrow(y)

  for (rep_i in seq_len(nrepeat)) {
    fold_id <- sample(rep(seq_len(folds), length.out = n))

    for (fold_i in seq_len(folds)) {
      test_idx <- which(fold_id == fold_i)
      train_idx <- setdiff(seq_len(n), test_idx)

      X_train <- lapply(X, function(b) b[train_idx, , drop = FALSE])
      X_test  <- lapply(X, function(b) b[test_idx, , drop = FALSE])

      y_train <- y[train_idx, , drop = FALSE]
      y_test  <- y[test_idx, , drop = FALSE]

      fit <- block.pls(
        X_train,
        y_train,
        ncomp = ncomp_max,
        mode = "regression"
      )

      pred <- predict(fit, newdata = X_test)
      pred_train <- predict(fit, newdata = X_train)
      # Usually pred$predict is an array: samples x yvars x comps
      pred_array <- pred$WeightedPredict
      pred_array_train <- pred_train$WeightedPredict

      for (comp_i in seq_len(ncomp_max)) {
        y_pred <- pred_array[, , comp_i, drop = FALSE]
        y_pred_train <- pred_array_train[, , comp_i, drop = FALSE]

        obs_test <- as.numeric(y_test[, 1])
        est_test <- as.numeric(y_pred[, 1, 1])
        obs_train <- as.numeric(y_train[, 1])
        est_train <- as.numeric(y_pred_train[, 1, 1])

        ok <- complete.cases(obs_test, est_test)

        if (sum(ok) < 3) {
          rmse_val <- NA_real_
          r2_val <- NA_real_
          q2_val <- NA_real_
        } else {
          obs_test_ok <- obs_test[ok]
          est_test_ok <- est_test[ok]

          rmse_val <- sqrt(mean((obs_test_ok - est_test_ok)^2))

          ss_res <- sum((obs_test_ok - est_test_ok)^2)
          ss_tot <- sum((obs_test_ok - mean(obs_train))^2)

          q2_val <- if (ss_tot == 0) NA_real_ else 1 - ss_res / ss_tot

          ss_res <- sum((obs_train - est_train)^2)
          ss_tot <- sum((obs_train - mean(obs_train))^2)

          r2_val <- if (ss_tot == 0) NA_real_ else 1 - ss_res / ss_tot


        }

        all_results[[res_i]] <- data.frame(
          repeat_id = rep_i,
          fold = fold_i,
          comp = comp_i,
          rmse = rmse_val,
          r2 = r2_val,
          q2 = q2_val
        )
        res_i <- res_i + 1

        print(paste0("comp ", comp_i, " of fold ", fold_i, " of rep ", rep_i, " done"))
      }
    }
  }

  all_results <- do.call(rbind, all_results)

  summary_res <- aggregate(
    cbind(rmse, r2, q2) ~ comp,
    data = all_results,
    FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
  )

  # flatten aggregate output
  summary_df <- data.frame(
    comp = summary_res$comp,
    rmse_mean = summary_res$rmse[, "mean"],
    rmse_sd   = summary_res$rmse[, "sd"],
    r2_mean = summary_res$r2[, "mean"],
    r2_sd   = summary_res$r2[, "sd"],
    q2_mean = summary_res$q2[, "mean"],
    q2_sd   = summary_res$q2[, "sd"]
  )

  best_comp_q2 <- summary_df$comp[which.max(summary_df$q2_mean)]
  best_comp_rmse <- summary_df$comp[which.min(summary_df$rmse_mean)]

  list(
    fold_results = all_results,
    summary = summary_df,
    best_comp_q2 = best_comp_q2,
    best_comp_rmse = best_comp_rmse
  )
}

run_modeling <- function (prot, met, img){
  met <- met[, colMeans(is.na(met)) <= NA_threshold]
  met <- met[rowMeans(is.na(met)) <= NA_threshold,]

  prot <- prot[, colMeans(is.na(prot)) <= NA_threshold]
  prot <- prot[rowMeans(is.na(prot)) <= NA_threshold,]

  img <- img[,colMeans(is.na(img)) <= NA_threshold]
  img <- img[rowMeans(is.na(img)) <= NA_threshold,]

  age_data <- na.omit(age_data)
  height_data <- na.omit(height_data)
  meta_data <- merge(age_data, height_data, by="eid")

  keep_eids <- Reduce(intersect, list(prot$eid, met$eid, img$eid, meta_data$eid))

  prot <- prot[match(keep_eids, prot$eid), ]
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

  prot_eids <- prot$eid

  meta_prot_aligned <- meta_data[match(prot$eid, meta_data$eid), ]

  prot <- prot[, setdiff(names(prot), "eid"), drop = FALSE]

  prot <- as.data.frame(lapply(prot, function(v) {
    residualize_vector(v, meta_prot_aligned, correction_formula)
  }))

  prot$eid <- prot_eids

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

  norm_prot <- fit_normalizer(prot, "eid")
  prot <- apply_normalizer(prot, norm_prot)

  norm_img <- fit_normalizer(img, "eid")
  img <- apply_normalizer(img, norm_img)

  X <- list(
    proteomics = prot[, setdiff(names(prot), "eid"), drop = FALSE],
    metabolomics = met[, setdiff(names(met), "eid"), drop = FALSE]
  )
  y <- img[[interested_feature]]

  univar_stats_protein <- data.frame(
    protein = colnames(X$proteomics),
    cor = NA_real_,
    pvalue = NA_real_
  )

  for (j in seq_along(X$proteomics)) {
    test_X <- X$proteomics[[j]]
    test_y <- y
    ok <- complete.cases(test_X, test_y)

    test <- suppressWarnings(cor.test(test_X[ok], test_y[ok], method = "pearson"))
    univar_stats_protein$cor[j] <- unname(test$estimate)
    univar_stats_protein$pvalue[j] <- test$p.value
  }

  univar_stats_met <- data.frame(
    metabolite = colnames(X$metabolomics),
    cor = NA_real_,
    pvalue = NA_real_
  )

  for (j in seq_along(X$metabolomics)) {
    test_X <- X$metabolomics[[j]]
    test_y <- y
    ok <- complete.cases(test_X, test_y)

    test <- suppressWarnings(cor.test(test_X[ok], test_y[ok], method = "pearson"))
    univar_stats_met$cor[j] <- unname(test$estimate)
    univar_stats_met$pvalue[j] <- test$p.value
  }

  prot[is.na(prot)] <- 0
  met[is.na(met)] <- 0
  img <- na.omit(img)

  keep_eids <- Reduce(intersect, list(prot$eid, met$eid, img$eid, meta_data$eid))

  met <- met[match(keep_eids, met$eid), ]
  prot <- prot[match(keep_eids, prot$eid), ]
  img <- img[match(keep_eids, img$eid), ]

  X <- list(
    proteomics = prot[, setdiff(names(prot), "eid"), drop = FALSE],
    metabolomics = met[, setdiff(names(met), "eid"), drop = FALSE]
  )
  row.names(X$proteomics) <- prot$eid
  row.names(X$metabolomics) <- met$eid

  y <- matrix(img[[interested_feature]], ncol=1)
  colnames(y) <- interested_feature
  rownames(y) <- keep_eids
  keep_proteins <- univar_stats_protein$protein[univar_stats_protein$pvalue < p_threshold]

  keep_met <- univar_stats_met$metabolite[univar_stats_met$pvalue < p_threshold]

  X_new <- list(
    proteomics = X$proteomics[, keep_proteins, drop = FALSE],
    metabolomics = X$metabolomics[, keep_met, drop = FALSE]
  )

  print(paste("Number of proteins selected:", ncol(X_new$proteomics)))
  print(paste("Number of metabolites selected:", ncol(X_new$metabolomics)))


  block.pls.result <- block.pls(X_new, y, ncomp=min(10, dim(X_new$proteomics), dim(X_new$metabolomics)), mode = 'regression', design="full")

  saveRDS(block.pls.result, paste0(interested_feature, "_results.rds"))

  png(filename=paste0(interested_feature, "_sample_spread.png"))
  plotIndiv(block.pls.result, ind.names = FALSE)
  dev.off()

  png(filename=paste0(interested_feature, "_Variance.png"))
  plotVar(block.pls.result, var.names = FALSE)
  dev.off()

  png(filename=paste0(interested_feature, "_Loadings_prot.png"))
  plotLoadings(block.pls.result, comp = 1, contrib = 'max', method = 'median', block="proteomics", title="Top 10 Proteins",ndisplay =10)
  dev.off()

  png(filename=paste0(interested_feature, "_Loadings_met.png"))
  plotLoadings(block.pls.result, comp = 1, contrib = 'max', method = 'median', block="metabolomics", title="Top 10 Proteins",ndisplay =10)
  dev.off()

  ncomp_max <- min(
    10,
    ncol(X_new$proteomics),
    ncol(X_new$metabolomics),
    nrow(y) - 1
  )

  cv_res <<- cv_block_pls_ncomp(
    X = X_new,
    y = y,
    folds = 10,
    nrepeat = 5,
    ncomp_max = ncomp_max,
    seed = 123
  )

  saveRDS(cv_res, paste0(interested_feature, "_R2_results.rds"))

  plot_R2 <- function(title, values, sd_values, ylab){
    vals <- values
    sd_vals <- sd_values
    best_idx <- which.max(vals)
    best_val <- vals[best_idx]
    x <- seq_along(vals)
    plot(vals,
         type="p",
         ylim = c(min(min(vals[vals > -1], 0)),
                       max(vals)),
         main=title,
        xlab="Components",
        ylab=ylab)

    grid()

    arrows(x, vals - sd_vals,
         x, vals + sd_vals,
         angle = 90, code = 3, length = 0.05)

    points(best_idx, best_val, col = "red", pch = 19, cex = 1.5)
    abline(h = best_val, col = "red", lty = 2)
    abline(h = 0, col = "red", lty = 1)
  }
  png(filename=paste0(interested_feature, "_R2_train.png"))
  plot_R2("R2_train_test", cv_res$summary$r2_mean, cv_res$summary$r2_sd, "R2")
  dev.off()

  png(filename=paste0(interested_feature, "_R2_heldout.png"))
  plot_R2("R2_heldout", cv_res$summary$q2_mean, cv_res$summary$q2_sd, "R2")
  dev.off()

  png(filename=paste0(interested_feature, "_R2_train_vs_test.png"))
    train_vals <- cv_res$summary$r2_mean
    train_sd   <- cv_res$summary$r2_sd

    test_vals  <- cv_res$summary$q2_mean   # or r2_test_mean if you have it
    test_sd    <- cv_res$summary$q2_sd

    x <- seq_along(train_vals)

    best_idx <- which.max(test_vals)
    best_val <- test_vals[best_idx]

    plot(x, train_vals,
         type = "p",
         ylim = c(min(c(train_vals - train_sd, test_vals - test_sd), na.rm = TRUE),
                  max(c(train_vals + train_sd, test_vals + test_sd), na.rm = TRUE)),
         main = "Train vs Heldout R2",
         xlab = "Components",
         ylab = "R2",
         col = "blue",
         xaxt = "n")

    axis(1, at = x)

    abline(v = x, col = "lightgray", lty = 3)
    abline(h = pretty(c(train_vals, test_vals), na.rm = TRUE),
           col = "lightgray", lty = 3)

    arrows(x, train_vals - train_sd,
           x, train_vals + train_sd,
           angle = 90, code = 3, length = 0.05, col = "blue")

    points(x, test_vals, col = "darkgreen")

    arrows(x, test_vals - test_sd,
           x, test_vals + test_sd,
           angle = 90, code = 3, length = 0.05, col = "darkgreen")

    lines(x, train_vals, col = "blue")
    lines(x, test_vals, col = "darkgreen")

    points(best_idx, best_val, col = "red", pch = 19, cex = 1.5)
    abline(v = best_idx, col = "red", lty = 3)
    abline(h = best_val, col = "red", lty = 2)

    legend("bottomright",
           legend = c("Train", "Heldout", paste0("Best component. R2=", round(best_val, 4))),
           col = c("blue", "darkgreen", "red"),
           pch = c(16, 16, 16),
           lty = c(1, 1, 2))

  dev.off()
}

run_modeling(proteomics_train, metabolomics_train, image_data_train)