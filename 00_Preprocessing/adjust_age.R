age_data <- read.csv("age_data.csv")


age_data <- na.omit(age_data)

age_data$age_between <- age_data$age_2 - age_data$age_0

write.csv(age_data, "age_data.csv", row.names = F)