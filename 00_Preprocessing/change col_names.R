data <- read.csv("age_data.csv")


colnames(data) <- c("eid", "age_0", "age_2")

write.csv(data, "age_data.csv", row.names = F)