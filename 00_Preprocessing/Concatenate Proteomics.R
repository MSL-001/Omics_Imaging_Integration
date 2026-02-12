prot1 <- read.csv("prot1.csv")
prot2 <- read.csv("prot2.csv")
prot3 <- read.csv("prot3.csv")
prot4 <- read.csv("prot4.csv")
prot5 <- read.csv("prot5.csv")
prot6 <- read.csv("prot6.csv")


prot12 <- merge(prot1, prot2, by="eid")
prot123 <- merge(prot12, prot3, by="eid")
prot1234 <- merge(prot123, prot4, by="eid")
prot12345 <- merge(prot1234, prot5, by="eid")
prot_all <- merge(prot12345, prot6, by="eid")


write.csv(prot_all, "prot_all_test.csv", row.names=F)
