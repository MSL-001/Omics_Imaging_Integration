metabolomics <- read.csv("metabolomics_data.csv")
image_features <- read.csv("volume_data.csv")
age_data <- read.csv("age_data.csv")

met_eids <- metabolomics$eid
img_eids <- image_features$eid
keep_eids <- intersect(met_eids, img_eids)

metabolomics <- metabolomics[metabolomics$eid %in% keep_eids, ]
image_features <- image_features[image_features$eid %in% keep_eids, ]
age_data <- age_data[age_data$eid %in% keep_eids, ]

print(paste("Metabolomics: ", dim(metabolomics), sep=" "))

print(paste("Image features: ", dim(image_features), sep=" "))

metabolomics <- metabolomics[rowMeans(is.na(metabolomics)) <= 0.30,]

metabolomics <- metabolomics[, colMeans(is.na(metabolomics)) <= 0.30]

image_features <- image_features[rowMeans(is.na(image_features)) <= 0.30,]

image_features <- image_features[,colMeans(is.na(image_features)) <= 0.30]

normalize <- function(df, exclude){
  data_cols <- setdiff(colnames(df), exclude)
  df[data_cols] <- lapply(df[data_cols], function(x) {
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(x)
    (x - mean(x, na.rm = TRUE)) / s
  })
}

met_eids <- metabolomics$eid
img_eids <- image_features$eid
keep_eids <- intersect(met_eids, img_eids)

metabolomics <- metabolomics[metabolomics$eid %in% keep_eids, ]
image_features <- image_features[image_features$eid %in% keep_eids, ]
age_data <- age_data[age_data$eid %in% keep_eids, ]

print(paste("Metabolomics after > 30 NA remove: ", dim(metabolomics), sep=" "))

print(paste("Image features > 30 NA remove: ", dim(image_features), sep=" "))

print(paste("Meta data > 30 NA remove: ", dim(age_data), sep=" "))

metabolomics <- na.omit(metabolomics)

image_features <- na.omit(image_features)

met_eids <- metabolomics$eid
img_eids <- image_features$eid
keep_eids <- intersect(met_eids, img_eids)

metabolomics <- metabolomics[metabolomics$eid %in% keep_eids, ]
image_features <- image_features[image_features$eid %in% keep_eids, ]
age_data <- age_data[age_data$eid %in% keep_eids, ]

print(paste("Metabolomics after full NA-remove: ", dim(metabolomics), sep=" "))

print(paste("Image features after full NA-remove: ", dim(image_features), sep=" "))


normalize(metabolomics, "eid")
normalize(image_features, c("eid", "bmi", "sex"))

for (feature in colnames(image_features)){
  if (!(feature %in% c("eid", "bmi", "sex"))){
    y <- image_features[[feature]]
    model <- lm(y ~ ., data = df)
    summary(model)
  }
}

