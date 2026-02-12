path_female_img <- "features_UKBB_vibe_80/ids_95k_female.txt"

path_male_img <- "features_UKBB_vibe_80/ids_95k_male.txt"

path_prot <- "Data/00_Preprocessing/proteomics_ids.csv"

path_met <- "Data/00_Preprocessing/metabolomics_ids.csv"




female_img_participants <- 
  read.table(path_female_img, sep=" ",header=F)

colnames(female_img_participants) <- "ids"

male_img_participants <- 
  read.table(path_male_img,sep=" ",header=F)

colnames(male_img_participants) <- "ids"

prot_participants <-
  read.table(path_prot, sep=" ",header=T)

colnames(prot_participants) <- "ids"

met_participants <-
  read.table(path_met, sep=" ",header=T)

colnames(met_participants) <- "ids"

met_img_male <- intersect(met_participants$ids, male_img_participants$ids)
met_img_female <- intersect(met_participants$ids, female_img_participants$ids)

prot_img_male <- intersect(prot_participants$ids, male_img_participants$ids)
prot_img_female <- intersect(prot_participants$ids, female_img_participants$ids)

prot_met_img_male <- intersect(prot_img_male, met_img_male)
prot_met_img_female <- intersect(prot_img_female, met_img_female)

sample_prot_male <- data.frame(sample(prot_img_male, 1000, replace=F))
sample_prot_female <- data.frame(sample(prot_img_female, 1000, replace=F))
sample_met_male <- data.frame(sample(prot_img_male, 1000, replace=F))
sample_met_female <- data.frame(sample(prot_img_female, 1000, replace=F))

colnames(sample_prot_male) <- "eid"
colnames(sample_prot_female) <- "eid"
colnames(sample_met_male) <- "eid"
colnames(sample_met_female) <- "eid"


write.csv(sample_prot_male, "Data/00_Preprocessing/pilot2_m_prot_cohort.csv", row.names=F)
write.csv(sample_prot_female, "Data/00_Preprocessing/pilot2_f_prot_cohort.csv", row.names=F)
write.csv(sample_met_male, "Data/00_Preprocessing/pilot2_m_met_cohort.csv", row.names=F)
write.csv(sample_met_female, "Data/00_Preprocessing/pilot2_f_met_cohort.csv", row.names=F)

