if(!require("ppcor")){
  install.packages("ppcor", repos='http://cran.us.r-project.org')
  library(ppcor)
}

calc_statistics <- function(img1, img2, df, covars) {
  cols_needed <- c(img1, img2, covars)
  d <- df[, cols_needed, drop = FALSE]
  ok <- complete.cases(d)

  n_used <- sum(ok)
  if (n_used > 20) {
    corr <- cor(df[[img1]], df[[img2]], use = "complete.obs")

    d <- d[ok, , drop = FALSE]

    correlations <- pcor(d)

    adj_corr <- correlations$estimate[img1, img2]

    p <- correlations$p.value[img1, img2]

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

image_features <- read.csv("features_UKBB_vibe_80/sex_bmi_meanvol.csv")
colnames(image_features)[colnames(image_features) == "sbj_id"] <- "eid"

img_cols  <- setdiff(colnames(image_features), c("eid", "sex", "bmi"))

joined <- image_features

joined[img_cols] <- lapply(joined[img_cols], function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x)
  (x - mean(x, na.rm = TRUE)) / s
})





img1_feat <- c()
img2_feat <- c()
corr_val <- c()
samples_used <- c()
covariate_adjusted <- c()
p_value <- c()



for (img1 in img_cols){
  for (img2 in img_cols){
    if (img1 != img2){
      res <- calc_statistics(img1, img2, joined)
      if (!is.na(res$adj_corr)){
        img1_feat <- c(img1_feat, img1)
        img2_feat <- c(img2_feat, img2)
        corr_val <- c(corr_val, res$corr)
        covariate_adjusted <- c(covariate_adjusted, res$adj_corr)
        samples_used <- c(samples_used, res$n)
        p_value <- c(p_value, res$p)
      }

    }
  }

}


df <- data.frame(
  "BC_Feature_1" = img1_feat,
  "BC_feature_2" = img2_feat,
  "Raw_corr" = corr_val,
  "Adjusted_corr" = covariate_adjusted,
  "P_value" = p_value,
  "Rows_used" = samples_used
)

df$q_value <- p.adjust(df$`P_value`, method = "BH")


write.csv(df, "correlations_adj.csv", row.names = F)

