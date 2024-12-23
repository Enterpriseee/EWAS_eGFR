#library soft ware
library(data.table)
library(dplyr)

#load data
a <- fread("data.csv")

#Calculate the mean and standard deviation of the control group
control_data <- a[a$corc == 0, c("CpG_9", "CpG_12","CpG_13","CpG_22")]
mean_control <- colMeans(control_data, na.rm = TRUE)
sd_control <- apply(control_data, 2, sd, na.rm = TRUE)

#Calculation of standardized values
standardized_data <- sweep(a[, c("CpG_9", "CpG_12","CpG_13","CpG_22")], 2, mean_control, "-")
standardized_data <- sweep(standardized_data, 2, sd_control, "/")

#Defining the direction of effect
W_c <- c(
  CpG_9 = -1, CpG_12 = -1, CpG_13 = -1, CpG_22 = -1
)
weighted_standardized_data <- sweep(standardized_data, 2, W_c, "*")

#Calculate MRS
MRS <- rowMeans(weighted_standardized_data, na.rm = TRUE)

#Adding results to the raw data
a$MRS <- MRS

#save data
write.csv(a,file = "data_with_MRS.csv",row.names = FALSE)
