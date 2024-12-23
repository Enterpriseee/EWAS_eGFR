rm(list = ls())
options(stringsAsFactors = FALSE)

#library software and load data
library(plyr)
library(doMC)
library(openxlsx)
doMC::registerDoMC(cores = 14)

phenfile <- 'sample_information.xlsx'
datafile <- 'rrbs_clean_data/predictedMeth_m.RDS'
refactor_obj_name <- 'rrbs_clean_data/refactor_obj.RDS'

phen <- read.xlsx(phenfile)
data_d <- readRDS(datafile)
data_d <- as.data.frame(data_d)

#Data processing: Ensure that phenotypic data correspond to methylation data
rownames(phen) <- phen$no
phen <- phen[colnames(data_d), ]
all(rownames(phen) == colnames(data_d))

#Preliminary GEE model fitting
age <- phen$age
gender <- as.factor(phen$gender)
fid <- as.factor(phen$family_ID)
eGFR <- phen$eGFR
PCs <- readRDS(refactor_obj_name)

require(gdata)
require(geepack)

dummy <- as.numeric(data_d[1, ])
dummy <- summary(geeglm(dummy ~ eGFR + age + gender + PCs$PC1 + PCs$PC2 + PCs$PC3 + PCs$PC4 + PCs$PC5, id = fid, corstr = "exchangeable"))$coefficient
dummy <- unmatrix(dummy, byrow = TRUE)
dummy <- cbind.data.frame(unit = rownames(data_d[1, ]), t(dummy))
dummy[1, ] <- NA

#Subject analysis: Parallel processing GEE model
require(plyr)
require(doMC)
doMC::registerDoMC(cores = 14)
require(gdata)
require(data.table)
result <- rbindlist(alply(data_d, 1, function(obs) {
  tryCatch({
    sumt <- summary(geeglm(as.numeric(obs) ~ eGFR + age + gender + PCs$PC1 + PCs$PC2 + PCs$PC3 + PCs$PC4 + PCs$PC5, id = fid, corstr = "exchangeable"))$coefficient
    sumt <- unmatrix(sumt, byrow = TRUE)
    sumt <- cbind.data.frame(unit = rownames(obs), t(sumt))
    sumt$unit <- rownames(obs)
    return(sumt)
  }, error = function(e) {
    dummyresult <- dummy
    dummyresult$unit <- rownames(obs)
    return(dummyresult)
  })
}, .progress = "none", .parallel = TRUE))

#save data
result <- as.data.frame(result)
saveRDS(result, 'rrbs_clean_data/lmer.RDS')

