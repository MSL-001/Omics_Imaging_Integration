args <- commandArgs(trailingOnly = TRUE)

input_omic_path <- args[1]
input_image_path <- args[2]
output_name <- args[3]
cohort_eids_path <- args[4]
gender <- args[5]


if(!require("ppcor")){
  install.packages("ppcor", repos='http://cran.us.r-project.org')
  library(ppcor)
}

library(dplyr)

calc_statistics <- function(omic, img, df, covars) {
  cols_needed <- c(omic, img, covars)
  d <- df[, cols_needed, drop = FALSE]
  ok <- complete.cases(d)

  n_used <- sum(ok)
  if (n_used >= 20) {
    corr_raw_pearson <- cor(df[[omic]], df[[img]], use = "complete.obs")
    corr_raw_spearman <- cor(df[[omic]], df[[img]], use = "complete.obs", method = "spearman")


    d <- d[ok, , drop = FALSE]

    pt_pearson <- tryCatch(
      ppcor::pcor.test(d[[omic]], d[[img]], d[, covars, drop = FALSE]),
      error = function(e) {
        message("FAIL: ", omic, " vs ", img, "  n=", n_used,
                "  sd(omic)=", sd(d[[omic]]), " sd(img)=", sd(d[[img]]),
                "  sd(covars)=", paste(round(sapply(d[, covars, drop = FALSE], sd), 6), collapse=", "),
                "  err=", conditionMessage(e))
      }
    )

    pt_spearman <- tryCatch(
      ppcor::pcor.test(d[[omic]], d[[img]], d[, covars, drop = FALSE], method="spearman"),
      error = function(e) {
        message("FAIL: ", omic, " vs ", img, "  n=", n_used,
                "  sd(omic)=", sd(d[[omic]]), " sd(img)=", sd(d[[img]]),
                "  sd(covars)=", paste(round(sapply(d[, covars, drop = FALSE], sd), 6), collapse=", "),
                "  err=", conditionMessage(e))
      }
    )


    return(list(corr_raw_pearson=corr_raw_pearson,
                corr_adj_pearson = if (is.null(pt_pearson)) NA else unname(pt_pearson$estimate),
                p_pearson = if (is.null(pt_pearson)) NA else unname(pt_pearson$p.value),
                n=n_used,
                corr_raw_spearman = corr_raw_spearman,
                corr_adj_spearman = if (is.null(pt_spearman)) NA else unname(pt_spearman$estimate),
                p_spearman = if (is.null(pt_spearman)) NA else unname(pt_spearman$p.value)))

  } else {
   return(list(corr_raw_pearson=NA,
                corr_adj_pearson=NA,
                p_pearson=NA,
                n=n_used,
                corr_raw_spearman = NA,
                corr_adj_spearman=NA,
                p_spearman=NA))
  }

}

eids <- read.csv(cohort_eids_path)
omics <- read.csv(input_omic_path, na.strings = "NA")
omics <- omics %>%
  filter(eid %in% eids$eid)

img_features <- read.csv(input_image_path, na.strings = "NA")
img_features <- img_features %>%
  filter(subject_ID %in% eids$eid)
age_data <- read.csv("age_data.csv")
height_data <- read.csv("height_data.csv")

joinedx <- merge(age_data, omics, by="eid")
joinedx <- merge(height_data, joinedx, by="eid")
joined <- merge(joinedx, img_features, by.x="eid", by.y="subject_ID")

omic_cols <- setdiff(colnames(omics), "eid")
img_cols  <- setdiff(colnames(img_features), "subject_ID")

covars_1 <- c("age_0", "age_between")

covars_2 <- c("age_0", "age_between", "p50_i0")

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

joined[covars_2] <- lapply(joined[covars_2], function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x)
  (x - mean(x, na.rm = TRUE)) / s
})

n_pairs <- length(omic_cols) * length(img_cols)

omics_feat <- character(n_pairs)
img_feat   <- character(n_pairs)
corr_pearson   <- numeric(n_pairs)
adj_corr_pearson    <- numeric(n_pairs)
p_pearson      <- numeric(n_pairs)
n_used     <- integer(n_pairs)
corr_spearman   <- numeric(n_pairs)
adj_corr_spearman    <- numeric(n_pairs)
p_spearman      <- numeric(n_pairs)
corr_pearson_height   <- numeric(n_pairs)
adj_corr_pearson_height    <- numeric(n_pairs)
p_pearson_height      <- numeric(n_pairs)
corr_spearman_height   <- numeric(n_pairs)
adj_corr_spearman_height    <- numeric(n_pairs)
p_spearman_height      <- numeric(n_pairs)

k <- 0L
for (omic in omic_cols){
  for (feature in img_cols){
    k <- k + 1L
    res <- calc_statistics(omic, feature, joined, covars_1)
    omics_feat[k] <- omic
    img_feat[k]   <- feature
    corr_pearson[k]   <- res$corr_raw_pearson
    adj_corr_pearson[k]    <- res$corr_adj_pearson
    p_pearson[k]      <- res$p_pearson
    n_used[k]     <- res$n
    corr_spearman[k] <- res$corr_raw_spearman
    adj_corr_spearman[k]    <- res$corr_adj_spearman
    p_spearman[k]      <- res$p_spearman
  }
}

k <- 0L
for (omic in omic_cols){
  for (feature in img_cols){
    k <- k + 1L
    res <- calc_statistics(omic, feature, joined, covars_2)
    corr_pearson_height[k]   <- res$corr_raw_pearson
    adj_corr_pearson_height[k]    <- res$corr_adj_pearson
    p_pearson_height[k]      <- res$p_pearson
    corr_spearman_height[k] <- res$corr_raw_spearman
    adj_corr_spearman_height[k]    <- res$corr_adj_spearman
    p_spearman_height[k]      <- res$p_spearman
  }
}

keep <- !is.na(adj_corr_pearson)

df <- data.frame(
  Omics_Feature  = omics_feat[keep],
  BC_feature     = img_feat[keep],
  Rows_used      = n_used[keep],
  Raw_corr_pearson       = corr_pearson[keep],
  Adjusted_corr_pearson  = adj_corr_pearson[keep],
  P_value_pearson        = p_pearson[keep],
  Raw_corr_spearman       = corr_spearman[keep],
  Adjusted_corr_spearman  = adj_corr_spearman[keep],
  P_value_spearman        = p_spearman[keep],
  Raw_corr_pearson_height       = corr_pearson_height[keep],
  Adjusted_corr_pearson_height  = adj_corr_pearson_height[keep],
  P_value_pearson_height        = p_pearson_height[keep],
  Raw_corr_spearman_height       = corr_spearman_height[keep],
  Adjusted_corr_spearman_height  = adj_corr_spearman_height[keep],
  P_value_spearman_height        = p_spearman_height[keep],
  stringsAsFactors = FALSE
)

df$q_value_pearson <- p.adjust(df$`P_value_pearson`, method = "BH")
df$q_value_spearman <- p.adjust(df$`P_value_spearman`, method = "BH")
df$q_value_pearson_height <- p.adjust(df$`P_value_pearson_height`, method = "BH")
df$q_value_spearman_height <- p.adjust(df$`P_value_spearman_height`, method = "BH")


# df <- df[, names(df) != "P_value", drop = FALSE]

write.csv(df, output_name, row.names = F)

