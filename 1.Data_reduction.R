#library software
library(data.table)
library(openxlsx)
library(dplyr)
library(tibble)

#load data
a<-read.xlsx("data.xlsx")

#prosess data
a <- a %>%
mutate_all(~ifelse(. == 0.00, 0.01, .))
a <- a %>%
mutate_all(~ifelse(. == 1.00, 0.99, .))
a <- a %>%
column_to_rownames(var = names(a)[1])
a <- log2(a / (1 - a))
a<-t(a)
a<-as.data.frame(a)
a$no <- rownames(a)
a <- a[, c("no", setdiff(names(a), "no"))]
rownames(a)<-NULL
colnames(a) <- gsub("SYNGR3-10_", "", colnames(a))

#Add number
b<-fread("Matching result.csv")
b$no<-as.numeric(b$no)
a$no<-as.numeric(a$no)
b<-merge(b,a,by="no")

#Add basic information
a<-read.xlsx("basic information.xlsx")
b<-merge(b,a,by="no")

#data 
b <- subset(b, select = -c(CpG_10.11, CpG_25.26.27.28.29.30))
b <- as.data.frame(b)
b$CpG_7 <- b$CpG_6.7
names(b)[names(b) == "CpG_6.7"] <- "CpG_6"
b$CpG_13 <- b$CpG_12.13
names(b)[names(b) == "CpG_12.13"] <- "CpG_12"
b$CpG_20 <- b$CpG_19.20.21
b$CpG_21 <- b$CpG_19.20.21
names(b)[names(b) == "CpG_19.20.21"] <- "CpG_19"
b <- b[, c(1:9, 30, 10:29, 31:ncol(b))]
b <- b[, c(1:13, 31, 14:30, 32:ncol(b))]
b <- b[, c(1:20, 32, 21:31, 33:ncol(b))]
b <- b[, c(1:21, 33, 22:32)]


#save data
write.csv(b,file = "data.csv",row.names = FALSE)
