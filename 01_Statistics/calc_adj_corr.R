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
    corr <- cor(df[[omic]], df[[img]], use = "complete.obs")

    d <- d[ok, , drop = FALSE]

    pt <- tryCatch(
      ppcor::pcor.test(d[[omic]], d[[img]], d[, covars, drop = FALSE]),
      error = function(e) {
        message("FAIL: ", omic, " vs ", img, "  n=", n_used,
                "  sd(omic)=", sd(d[[omic]]), " sd(img)=", sd(d[[img]]),
                "  sd(covars)=", paste(round(sapply(d[, covars, drop = FALSE], sd), 6), collapse=", "),
                "  err=", conditionMessage(e))
      }
    )

    if (is.null(pt)){
      return(list(corr=corr, adj_corr=NA, p=NA, n=n_used))
    }

    return(list(corr=corr,
                adj_corr=unname(pt$estimate),
                p=unname(pt$p.value),
                n=n_used))

  } else {
    return(list(corr = NA,
                adj_corr = NA,
                p    = NA,
                n    = n_used)
    )
  }

}

eids <- read.csv(cohort_eids_path)
omics <- read.csv(input_omic_path, na.strings = "NA")
omics <- omics %>%
  filter(eid %in% eids$eid)
print(dim(omics))

img_features <- read.csv(input_image_path, na.strings = "NA")
img_features <- img_features %>%
  filter(subject_ID %in% eids$eid)
age_data <- read.csv("age_data.csv")
img_features <- img_features %>%
  filter(subject_ID %in% eids$eid)
print(dim(img_features))


joinedx <- merge(age_data, omics, by="eid")
joined <- merge(joinedx, img_features, by.x="eid", by.y="subject_ID")

print(dim(joined))

omic_cols <- setdiff(colnames(omics), "eid")
img_cols  <- setdiff(colnames(img_features), "subject_ID")

covars <- c("age_0", "age_between")

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

joined[covars] <- lapply(joined[covars], function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x)
  (x - mean(x, na.rm = TRUE)) / s
})

n_pairs <- length(omic_cols) * length(img_cols)

omics_feat <- character(n_pairs)
img_feat   <- character(n_pairs)
corr_val   <- numeric(n_pairs)
adj_val    <- numeric(n_pairs)
p_val      <- numeric(n_pairs)
n_used     <- integer(n_pairs)

k <- 0L
for (omic in omic_cols){
  for (feature in img_cols){
    k <- k + 1L
    res <- calc_statistics(omic, feature, joined, covars)
    omics_feat[k] <- omic
    img_feat[k]   <- feature
    corr_val[k]   <- res$corr
    adj_val[k]    <- res$adj_corr
    p_val[k]      <- res$p
    n_used[k]     <- res$n
  }
}

keep <- !is.na(adj_val)

df <- data.frame(
  Omics_Feature  = omics_feat[keep],
  BC_feature     = img_feat[keep],
  Raw_corr       = corr_val[keep],
  Adjusted_corr  = adj_val[keep],
  P_value        = p_val[keep],
  Rows_used      = n_used[keep],
  stringsAsFactors = FALSE
)

df$q_value <- p.adjust(df$`P_value`, method = "BH")

df <- df[, names(df) != "P_value", drop = FALSE]

write.csv(df, output_name, row.names = F)

