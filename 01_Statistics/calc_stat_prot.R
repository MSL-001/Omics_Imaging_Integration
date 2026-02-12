if(!require("ppcor")){
  install.packages("ppcor", repos='http://cran.us.r-project.org')
  library(ppcor)
}

calc_statistics <- function(omic, img, df, covars) {
  cols_needed <- c(omic, img, covars)
  d <- df[, cols_needed, drop = FALSE]
  ok <- complete.cases(d)

  n_used <- sum(ok)
  if (n_used > 20) {
    corr <- cor(df[[omic]], df[[img]], use = "complete.obs")

    d <- d[ok, , drop = FALSE]

    correlations <- pcor(d)

    adj_corr <- correlations$estimate[omic, img]

    p <- correlations$p.value[omic, img]

    return(list(corr = corr,
                adj_corr = adj_corr,
                p    = p,
                n    = n_used)
    )
  } else {
    return(list(corr = NA,
                adj_corr = NA,
                p    = NA,
                n    = n_used)
    )
  }

}

proteomics <- read.csv("proteomics_data.csv")
image_features <- read.csv("volume_data.csv")
age_data <- read.csv("age_data.csv")

joinedx <- merge(age_data, proteomics, by="eid")
joined <- merge(joinedx, image_features, by="eid")

omic_cols <- setdiff(colnames(proteomics), "eid")
img_cols  <- setdiff(colnames(image_features), c("eid", "sex", "bmi"))
covars <- c("age_0", "age_2")

joined[img_cols] <- lapply(joined[img_cols], function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x)
  (x - mean(x, na.rm = TRUE)) / s
})

joined[omic_cols] <- lapply(joined[omic_cols], function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x)
  (x - mean(x, na.rm = TRUE)) / s
})




omics_feat <- c()
img_feat <- c()
corr_val <- c()
samples_used <- c()
covariate_adjusted <- c()
p_value <- c()



for (omic in omic_cols){
  for (feature in img_cols){
    res <- calc_statistics(omic, feature, joined, covars)
    if (!is.na(res$adj_corr)){
      omics_feat <- c(omics_feat, omic)
      img_feat <- c(img_feat, feature)
      corr_val <- c(corr_val, res$corr)
      covariate_adjusted <- c(covariate_adjusted, res$adj_corr)
      samples_used <- c(samples_used, res$n)
      p_value <- c(p_value, res$p)
    }
  }

}


df <- data.frame(
  "Omics_Feature" = omics_feat,
  "BC_feature" = img_feat,
  "Raw_corr" = corr_val,
  "Adjusted_corr" = covariate_adjusted,
  "P_value" = p_value,
  "Rows_used" = samples_used
)

df$q_value <- p.adjust(df$`P_value`, method = "BH")


write.csv(df, "correlations_adj.csv", row.names = F)

