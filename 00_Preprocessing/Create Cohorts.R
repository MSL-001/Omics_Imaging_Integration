path_female_img <- "features_UKBB_vibe_80/ids_95k_female.txt"

path_male_img <- "features_UKBB_vibe_80/ids_95k_male.txt"

path_2prot <- "Data/00_Preprocessing/double_proteomics.csv"

path_1prot <- "Data/00_Preprocessing/proteomics_b0-6_eids.csv"

path_met <- "Data/00_Preprocessing/metabolomics_ids.csv"

set.seed(123)


female_img_participants <- 
  read.table(path_female_img, sep=" ",header=F)

colnames(female_img_participants) <- "eid"

male_img_participants <- 
  read.table(path_male_img,sep=" ",header=F)

colnames(male_img_participants) <- "eid"

double_prot_participants <-
  read.table(path_2prot, sep=" ",header=T)

colnames(double_prot_participants) <- "eid"

single_prot_participants <-
  read.table(path_1prot, sep=" ",header=T)

colnames(single_prot_participants) <- "eid"

met_participants <-
  read.table(path_met, sep=" ",header=T)

colnames(met_participants) <- "eid"

met_img_male <- intersect(met_participants$eid, male_img_participants$eid)
met_img_female <- intersect(met_participants$eid, female_img_participants$eid)

prot2_img_male <- intersect(double_prot_participants$eid, male_img_participants$eid)
prot2_img_female <- intersect(double_prot_participants$eid, female_img_participants$eid)

prot2_met_img_male <- intersect(prot2_img_male, met_img_male)
prot2_met_img_female <- intersect(prot2_img_female, met_img_female)

prot1_img_male <- intersect(single_prot_participants$eid, male_img_participants$eid)
prot1_img_female <- intersect(single_prot_participants$eid, female_img_participants$eid)

prot1_met_img_male <- intersect(prot1_img_male, met_img_male)
prot1_met_img_female <- intersect(prot1_img_female, met_img_female)


split_test_train <- function(eid_list, fraction){
  n <- length(eid_list)

  fold_idx <- sample(seq_len(n), size = floor(n*fraction), replace = FALSE)

  test <- data.frame(eid_list[fold_idx])
  train <- data.frame(eid_list[-fold_idx])
  return(list(test = test,
              train = train))
}

fraction <- 0.1

folds <- split_test_train(prot1_img_female, fraction)

prot1_img_female_test <- folds$test
prot1_img_female_train <- folds$train

folds <- split_test_train(prot1_img_male, fraction)

prot1_img_male_test <- folds$test
prot1_img_male_train <- folds$train

folds <- split_test_train(prot1_met_img_female, fraction)

prot1_met_img_female_test <- folds$test
prot1_met_img_female_train <- folds$train

folds <- split_test_train(prot1_met_img_male, fraction)

prot1_met_img_male_test <- folds$test
prot1_met_img_male_train <- folds$train

folds <- split_test_train(met_img_female, fraction)

met_img_female_test <- folds$test
met_img_female_train <- folds$train

folds <- split_test_train(met_img_male, fraction)

met_img_male_test <- folds$test
met_img_male_train <- folds$train

colnames(prot1_img_female_test) <- "eid"
colnames(prot1_img_female_train) <- "eid"
colnames(prot1_img_male_test) <- "eid"
colnames(prot1_img_male_train) <- "eid"
colnames(met_img_female_test) <- "eid"
colnames(met_img_female_train) <- "eid"
colnames(met_img_male_test) <- "eid"
colnames(met_img_male_train) <- "eid"
colnames(prot1_met_img_female_test) <- "eid"
colnames(prot1_met_img_female_train) <- "eid"
colnames(prot1_met_img_male_test) <- "eid"
colnames(prot1_met_img_male_train) <- "eid"


write.csv(prot1_img_female_test, "Data/00_Preprocessing/Single_prot_f/test_cohort.csv", row.names=F)
write.csv(prot1_img_female_train, "Data/00_Preprocessing/Single_prot_f/train_cohort.csv", row.names=F)
write.csv(prot1_img_male_test, "Data/00_Preprocessing/Single_prot_m/test_cohort.csv", row.names=F)
write.csv(prot1_img_male_train, "Data/00_Preprocessing/Single_prot_m/train_cohort.csv", row.names=F)

write.csv(prot1_met_img_female_test, "Data/00_Preprocessing/Met_prot_f/test_cohort.csv", row.names=F)
write.csv(prot1_met_img_female_train, "Data/00_Preprocessing/Met_prot_f/train_cohort.csv", row.names=F)
write.csv(prot1_met_img_male_test, "Data/00_Preprocessing/Met_prot_m/test_cohort.csv", row.names=F)
write.csv(prot1_met_img_male_train, "Data/00_Preprocessing/Met_prot_m/train_cohort.csv", row.names=F)

write.csv(met_img_female_test, "Data/00_Preprocessing/Met_f/test_cohort.csv", row.names=F)
write.csv(met_img_female_train, "Data/00_Preprocessing/Met_f/train_cohort.csv", row.names=F)
write.csv(met_img_male_test, "Data/00_Preprocessing/Met_m/test_cohort.csv", row.names=F)
write.csv(met_img_male_train, "Data/00_Preprocessing/Met_m/train_cohort.csv", row.names=F)
