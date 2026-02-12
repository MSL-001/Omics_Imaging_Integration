data <- read.csv("data.csv")


colnames(data) <- c("eid", "age_0", "age_2")

write.csv(data, "data_new.csv", row.names = F)