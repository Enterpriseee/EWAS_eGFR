rm(list = ls()) 

options(stringsAsFactors = FALSE)

#library data and load data
library(openxlsx)
sigresult <- read.xlsx ("sigresult.xlsx")

#extracting the CpG sites data
lmer_top_CpG <- subset(sigresult, pvalue < 0.00001)

#extracting the CpG sites number
myvars <- c("id_row_number")
lmer_top_CpG_unit <- lmer_top_CpG[myvars]

#rename
library(plyr)
lmer_top_CpG_unit <- rename(lmer_top_CpG_unit, c (id_row_number = "unit"))

# reading the methylation data of all CpG sites
datafile <- 'predictedMeth_m.RDS'
data_d <- readRDS(datafile)
data_d <- as.data.frame(data_d)
# adding the CpG sites name-unit
data_d$unit <- rownames(data_d)

# combine the data
library(plyr)
dat_CpG <- join(data_d, lmer_top_CpG_unit, by = 'unit', type = "right")

#rename the rownames
dat_CpG <- dat_CpG[!duplicated(dat_CpG$unit), ]
myvars1 <- dat_CpG[ , 1]
rownames(dat_CpG) <- myvars1

#delete the first column
myvars2 <- names(dat_CpG) %in% c("unit")
dat_CpG0 <- dat_CpG [!myvars2]

#save data
write.xlsx(dat_CpG0, "top_CpG_level.xlsx", rowNames = FALSE)
saveRDS(dat_CpG0,  "top_CpG_level.RDS")

