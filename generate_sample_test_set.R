set.seed(123)

# Generate Proteomics data
n_samples  <- 10000
n_proteins <- 500

proteomics_data <- data.frame(
  eid = 1001:(n_samples+1000)
)

proteomics_data[paste0("prot_", 1:n_proteins)] <-
  matrix(rnorm(n_samples * n_proteins, mean = 0, sd = 1),
         nrow = n_samples)


missing_rate <- 0.1

mask <- matrix(
  runif(n_samples * n_proteins) < missing_rate,
  nrow = n_samples
)

proteomics_data[-1][mask] <- NA

# Generate Metabolomics data
n_samples  <- 10000
n_metabolites <- 200

metabolomics_data <- data.frame(
  eid = 1001:(n_samples+1000)
)

metabolomics_data[paste0("met_", 1:n_metabolites)] <-
  matrix(rnorm(n_samples * n_metabolites, mean = 0, sd = 1),
         nrow = n_samples)


# missing_rate <- 0.1
#
# mask <- matrix(
#   runif(n_samples * n_metabolites) < missing_rate,
#   nrow = n_samples
# )
#
# metabolomics_data[-1][mask] <- NA

# Generate image data
n_samples <- nrow(proteomics_data)
n_features <- 20

generate_mock_meanff <- function(eids, n_features=20){

  n <- length(eids)

  df <- data.frame(subject_ID = eids)

  df[paste0("meanff_",1:n_features)] <-
    matrix(
      pmin(pmax(rnorm(n*n_features, 300, 120),0),1000),
      nrow=n
    )

  return(df)
}

image_data <- generate_mock_meanff(proteomics_data$eid)


train_idx <- sample(seq_len(n_samples), size = floor(n_samples * 0.9))

train_eids <- data.frame(proteomics_data$eid[train_idx])
test_eids  <- data.frame(proteomics_data$eid[-train_idx])

colnames(train_eids) <- "eid"
colnames(test_eids) <- "eid"

n <- length(train_eids)

fold_idx <- sample(seq_len(n), size = floor(n * 0.5))

train_eids_f1 <- train_eids[fold_idx, drop = FALSE]
train_eids_f2 <- train_eids[-fold_idx, drop = FALSE]

n_samples <- nrow(proteomics_data)

age_data <- data.frame(
  eid = proteomics_data$eid,
  age_0 = sample(45:70, n_samples, replace = TRUE),
  age_between = sample(6:14, n_samples, replace = TRUE)
)

height_data <- data.frame(
  eid = proteomics_data$eid,
  p50_i0 = sample(140:200, n_samples, replace = TRUE)
)


image_data$meanff_1 <- image_data$meanff_1 +
                        proteomics_data$prot_3 * 500 - proteomics_data$prot_1 *200 + proteomics_data$prot_5*100 + age_data$age_0 *1000

image_data$meanff_1 <- proteomics_data$prot_3 * 500 - proteomics_data$prot_1 *200 + proteomics_data$prot_5*100 + metabolomics_data$met_1*500
# image_data$meanff_1 <- proteomics_data$prot_3 * 500 - proteomics_data$prot_1 *200 + proteomics_data$prot_5*100



write.csv(train_eids, "train_cohort.csv", row.names=F)
write.csv(test_eids, "test_cohort.csv", row.names=F)