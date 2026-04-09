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
encoding_path <- "Data/Metadata/encoding.csv"

args <- commandArgs(trailingOnly = TRUE)

input_proteomics_path <- args[1]
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
proteomics_data <- read.csv(input_proteomics_path)
image_data <- read.csv(input_image_path)
age_data <- read.csv("age_data.csv")
encoding_path <- "encoding.csv"
height_data <- read.csv("height_data.csv")

encoding <- read.csv(encoding_path)




NA_threshold <- 0.3

# lookup <- setNames(encoding$Name, encoding$Code)
#
# colnames(proteomics_data) <- ifelse(
#   colnames(proteomics_data) %in% names(lookup),
#   lookup[colnames(proteomics_data)],
#   colnames(proteomics_data)
# )

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

fit_residualizer <- function(v, covar_data, correction_formula) {
  df <- covar_data
  df$target <- v

  fit <- lm(
    update(correction_formula, target ~ .),
    data = df,
    na.action = na.exclude
  )

  list(
    fit = fit,
    residuals = residuals(fit)
  )
}

apply_residualizer <- function(v, covar_data, fitted_obj) {
  df <- covar_data
  df$target <- v

  pred <- predict(fitted_obj$fit, newdata = df)
  v - pred
}

cv_block_pls_ncomp <- function(X, y, folds = 10, nrepeat = 5, ncomp_max = NULL, seed = 123) {
  set.seed(seed)

  if (is.vector(y) || is.null(dim(y))) {
    y <- data.frame(y = y)
  }

  all_results <- list()
  res_i <- 1

  n <- nrow(y)

  for (rep_i in seq_len(nrepeat)) {
    fold_id <- sample(rep(seq_len(folds), length.out = n))

    for (fold_i in seq_len(folds)) {
      test_idx <- which(fold_id == fold_i)
      train_idx <- setdiff(seq_len(n), test_idx)

      X_train <- X[train_idx, , drop = FALSE]
      X_test  <- X[test_idx, , drop = FALSE]

      y_train <- y[train_idx, , drop = FALSE]
      y_test  <- y[test_idx, , drop = FALSE]

      fit <- pls(
        X_train,
        y_train,
        ncomp = ncomp_max,
        mode = "regression"
      )

      pred <- predict(fit, newdata = X_test)
      pred_train <- predict(fit, newdata = X_train)
      pred_array <- pred$predict
      pred_array_train <- pred_train$predict

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

run_modeling <- function (prot, img){
  prot <- prot[, colMeans(is.na(prot)) <= NA_threshold]
  prot <- prot[rowMeans(is.na(prot)) <= NA_threshold,]

  img <- img[,colMeans(is.na(img)) <= NA_threshold]
  img <- img[rowMeans(is.na(img)) <= NA_threshold,]

  age_data <- na.omit(age_data)
  height_data <- na.omit(height_data)
  meta_data <- merge(age_data, height_data, by="eid")

  keep_eids <- Reduce(intersect, list(prot$eid, img$eid, meta_data$eid))

  prot <- prot[match(keep_eids, prot$eid), ]
  img <- img[match(keep_eids, img$eid), ]
  meta_data <- meta_data[match(keep_eids, meta_data$eid), ]


  prot_eids <- prot$eid

  meta_prot_aligned <- meta_data[match(prot$eid, meta_data$eid), ]

  prot <- prot[, setdiff(names(prot), "eid"), drop = FALSE]

  prot_train_mat <- prot[, setdiff(names(prot), "eid"), drop = FALSE]

  prot_resid_models <- vector("list", ncol(prot_train_mat))
  names(prot_resid_models) <- colnames(prot_train_mat)

  prot_resid <- prot_train_mat

  for (j in seq_along(prot_train_mat)) {
    res_obj <- fit_residualizer(
      v = prot_train_mat[[j]],
      covar_data = meta_prot_aligned,
      correction_formula = correction_formula
    )

    prot_resid[[j]] <- res_obj$residuals
    prot_resid_models[[j]] <- res_obj
  }

  prot_resid$eid <- prot_eids

  prot <- prot_resid

  img_eids <- img$eid

  img <- img[, interested_feature, drop = FALSE]

  meta_img_aligned <- meta_data[match(img_eids, meta_data$eid), ]

  img_train_fit <- fit_residualizer(
    v = img[[interested_feature]],
    covar_data = meta_img_aligned,
    correction_formula = correction_formula
  )

  img[[interested_feature]] <- img_train_fit$residuals

  img[["eid"]] <- img_eids

  saveRDS(prot_resid_models, paste0(interested_feature, "_protein_residualizers.rds"))
  saveRDS(img_train_fit,     paste0(interested_feature, "_image_residualizer.rds"))

  norm_prot <- fit_normalizer(prot, "eid")
  prot <- apply_normalizer(prot, norm_prot)

  norm_img <- fit_normalizer(img, "eid")
  img <- apply_normalizer(img, norm_img)

  saveRDS(norm_prot, paste0(interested_feature, "_protein_normalizers.rds"))
  saveRDS(norm_img, paste0(interested_feature, "_image_normalizer.rds"))


  X <- prot[, setdiff(names(prot), "eid"), drop = FALSE]
  y <- img[[interested_feature]]

  univar_stats <- data.frame(
    protein = colnames(X),
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

  prot[is.na(prot)] <- 0
  img <- na.omit(img)

  keep_eids <- Reduce(intersect, list(prot$eid, img$eid, meta_data$eid))

  prot <- prot[match(keep_eids, prot$eid), ]
  img <- img[match(keep_eids, img$eid), ]

  X <- prot[, setdiff(names(prot), "eid"), drop = FALSE]

  row.names(X) <- prot$eid

  y <- img[[interested_feature]]

  keep_proteins <- univar_stats$protein[univar_stats$pvalue < p_threshold]

  X <- X[, keep_proteins, drop = FALSE]

  print(paste("Number of proteins selected:", ncol(X)))

  pls.result <- pls(X, y, ncomp= 10, mode = 'regression')

  saveRDS(pls.result, paste0(interested_feature, "_results.rds"))

  png(filename=paste0(interested_feature, "_sample_spread.png"))
  plotIndiv(pls.result, ind.names = FALSE)
  dev.off()

  png(filename=paste0(interested_feature, "_Variance.png"))
  plotVar(pls.result, var.names = FALSE)
  dev.off()

  png(filename=paste0(interested_feature, "_Loadings.png"))
  plotLoadings(pls.result, comp = 1, contrib = 'max', method = 'median', block = "X", title="Top 10 Proteins",ndisplay =10)
  dev.off()

  Q2.pls.result <- perf(pls.result, validation = 'Mfold',
                folds = 10, nrepeat = 5)

  saveRDS(Q2.pls.result, paste0(interested_feature, "_R2_results.rds"))

  png(filename=paste0(interested_feature, "_R2_train_added_value.png"))
  vals <- Q2.pls.result$measures$R2$summary$mean
  p <- plot(Q2.pls.result, criterion = "R2", ylim=c(1,-1),title=paste0("R2 of comp 1 is ", round(Q2.pls.result$measures$R2$summary$mean[[1]], 4)))
  p <- p + coord_cartesian(ylim = c(min(min(vals[vals > -1], 0)), max(Q2.pls.result$measures$R2$summary$mean)))
  print(p)
  dev.off()

  png(filename=paste0(interested_feature, "_R2_heldout_added_value.png"))
  Q2.pls.result$measures$R2 <- Q2.pls.result$measures$Q2.total
  vals <- Q2.pls.result$measures$R2$summary$mean
  p <- plot(Q2.pls.result, criterion = "R2", ylim=c(1,-1),title=paste0("R2 of comp 1 is ", round(Q2.pls.result$measures$R2$summary$mean[[1]], 4)))
  p <- p + coord_cartesian(ylim = c(min(min(vals[vals > -1], 0)), max(Q2.pls.result$measures$R2$summary$mean)))
  print(p)
  dev.off()

  ncomp_max <- min(
    10,
    ncol(X),
    nrow(y) - 1
  )

  cv_res <- cv_block_pls_ncomp(
    X = X,
    y = y,
    folds = 10,
    nrepeat = 5,
    ncomp_max = ncomp_max,
    seed = 123
  )

  saveRDS(cv_res, paste0(interested_feature, "_R2_results_2.rds"))

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

run_modeling(proteomics_train, image_data_train)