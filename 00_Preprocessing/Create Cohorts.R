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

length(met_img_male)
length(met_img_female)

length(prot2_img_male)
length(prot2_img_female)
length(prot1_img_male)
length(prot1_img_female)

length(prot2_met_img_male)
length(prot2_met_img_female)
length(prot1_met_img_male)
length(prot1_met_img_female)

generate_folds <- function(eid_list){
  n <- length(eid_list)

  fold_idx <- sample(seq_len(n), size = floor(n/2), replace = FALSE)

  fold1 <- data.frame(eid_list[fold_idx])
  fold2 <- data.frame(eid_list[-fold_idx])
  return(list(fold1 = fold1,
              fold2 = fold2))
}
folds <- generate_folds(prot1_img_female)

prot1_img_female_fold1 <- folds$fold1
prot1_img_female_fold2 <- folds$fold2

folds <- generate_folds(prot1_img_male)

prot1_img_male_fold1 <- folds$fold1
prot1_img_male_fold2 <- folds$fold2

folds <- generate_folds(met_img_female)

met_img_female_fold1 <- folds$fold1
met_img_female_fold2 <- folds$fold2

folds <- generate_folds(met_img_male)

met_img_male_fold1 <- folds$fold1
met_img_male_fold2 <- folds$fold2

folds <- generate_folds(prot1_met_img_female)

prot1_met_img_female_fold1 <- folds$fold1
prot1_met_img_female_fold2 <- folds$fold2

folds <- generate_folds(prot1_met_img_male)

prot1_met_img_male_fold1 <- folds$fold1
prot1_met_img_male_fold2 <- folds$fold2

prot2_img_male <- data.frame(prot2_img_male)
prot2_img_female <- data.frame(prot2_img_female)


colnames(prot1_img_female_fold1) <- "eid"
colnames(prot1_img_female_fold2) <- "eid"
colnames(prot1_img_male_fold1) <- "eid"
colnames(prot1_img_male_fold2) <- "eid"
colnames(met_img_female_fold1) <- "eid"
colnames(met_img_female_fold2) <- "eid"
colnames(met_img_male_fold1) <- "eid"
colnames(met_img_male_fold2) <- "eid"
colnames(prot1_met_img_female_fold1) <- "eid"
colnames(prot1_met_img_female_fold2) <- "eid"
colnames(prot1_met_img_male_fold1) <- "eid"
colnames(prot1_met_img_male_fold2) <- "eid"
colnames(prot2_img_female) <- "eid"
colnames(prot2_img_male) <- "eid"

write.csv(prot1_img_female_fold1, "Data/00_Preprocessing/Single_prot_f/fold1/cohort.csv", row.names=F)
write.csv(prot1_img_female_fold2, "Data/00_Preprocessing/Single_prot_f/fold2/cohort.csv", row.names=F)

write.csv(prot1_img_male_fold1, "Data/00_Preprocessing/Single_prot_m/fold1/cohort.csv", row.names=F)
write.csv(prot1_img_male_fold2, "Data/00_Preprocessing/Single_prot_m/fold2/cohort.csv", row.names=F)

write.csv(prot2_img_female, "Data/00_Preprocessing/Double_prot_f/cohort.csv", row.names=F)
write.csv(prot2_img_male, "Data/00_Preprocessing/Double_prot_m/cohort.csv", row.names=F)


write.csv(met_img_female_fold1, "Data/00_Preprocessing/Met_f/fold1/cohort.csv", row.names=F)
write.csv(met_img_female_fold2, "Data/00_Preprocessing/Met_f/fold2/cohort.csv", row.names=F)

write.csv(met_img_male_fold1, "Data/00_Preprocessing/Met_m/fold1/cohort.csv", row.names=F)
write.csv(met_img_male_fold2, "Data/00_Preprocessing/Met_m/fold2/cohort.csv", row.names=F)

write.csv(prot1_met_img_female_fold1, "Data/00_Preprocessing/Met_prot_f/fold1/cohort.csv", row.names=F)
write.csv(prot1_met_img_female_fold2, "Data/00_Preprocessing/Met_prot_f/fold2/cohort.csv", row.names=F)

write.csv(prot1_met_img_male_fold1, "Data/00_Preprocessing/Met_prot_m/fold1/cohort.csv", row.names=F)
write.csv(prot1_met_img_male_fold2, "Data/00_Preprocessing/Met_prot_m/fold2/cohort.csv", row.names=F)
