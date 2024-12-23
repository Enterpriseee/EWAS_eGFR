options(stringsAsFactors = FALSE)

#library software and load data
library(openxlsx)
sample <- read.csv ("eGFRtraits.csv")

#create path
dir.create("ENSG")

#gee analysis
sample$age <- sample$age
sample$eGFR <- sample$eGFR
sample$gender <- as.factor(sample$gender)
sample$fid <- as.factor(sample$family_ID)
sample$eGFR <- as.numeric(sample$eGFR)

require(geepack)
fit <- geeglm(eGFR ~  age + gender + ENSG00000127561, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000127561.xlsx"))

require(geepack)
fit <- geeglm(eGFR ~  age + gender + ENSG00000196365, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000196365.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000115844, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000115844.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000185730, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000185730.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000188107, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000188107.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000185551, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000185551.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000104897, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000104897.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000139835, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000139835.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000100425, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000100425.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000165996, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000165996.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000130382, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000130382.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000189067, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000189067.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000175785, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000175785.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000180264, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000180264.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000148120, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000148120.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000108001, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000108001.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000182580, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000182580.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000066735, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000066735.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000164690, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000164690.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000184640, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000184640.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000160255, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000160255.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000171735, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000171735.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000164347, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000164347.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000204103, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000204103.xlsx"))

fit <- geeglm(eGFR ~  age + gender + ENSG00000165606, id = fid, corstr = "exchangeable", data = sample)
summary(fit)
summary_text <- capture.output(summary(fit))
writeLines(summary_text, file("ENSG/ENSG00000165606.xlsx"))

