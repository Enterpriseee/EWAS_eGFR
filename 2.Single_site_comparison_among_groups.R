#library software
library(epiDisplay)
library(survival)
library(Epi)
library(data.table)

#load data
data<-fread("data.csv")
data<-as.data.frame(data)

#Normal distribution test
results <- sapply(data[, 4:25], function(x) shapiro.test(x)$p.value)
results

#comparison among groups
#CpG1
res <- wilcox.test(CpG_1 ~ corc, data = data)
res
result <- data.frame(
  CpG = "CpG_1",    
  W = res$statistic,    
  P = res$p.value       
)
#CpG2
res <- wilcox.test(CpG_2 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_2",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG3
res <- wilcox.test(CpG_3 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_3",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG4
res <- wilcox.test(CpG_4 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_4",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG5
res <- wilcox.test(CpG_5 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_5",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG6
res <- wilcox.test(CpG_6 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_6",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG7
res <- wilcox.test(CpG_7 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_7",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG8
res <- wilcox.test(CpG_8 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_8",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG9
res <- wilcox.test(CpG_9 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_9",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG12
res <- wilcox.test(CpG_12 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_12",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG13
res <- wilcox.test(CpG_13 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_13",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG14
res <- wilcox.test(CpG_14 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_14",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG15
res <- wilcox.test(CpG_15 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_15",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG16
res <- wilcox.test(CpG_16 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_16",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG17
res <- wilcox.test(CpG_17 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_17",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG18
res <- wilcox.test(CpG_18 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_18",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG19
res <- wilcox.test(CpG_19 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_19",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG20
res <- wilcox.test(CpG_20 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_20",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG21
res <- wilcox.test(CpG_21 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_21",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG22
res <- wilcox.test(CpG_22 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_22",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG23
res <- wilcox.test(CpG_23 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_23",       
  W = res$statistic, 
  P = res$p.value    
))
#CpG24
res <- wilcox.test(CpG_24 ~ corc, data = data)
result <- rbind(result, data.frame(
  CpG = "CpG_24",       
  W = res$statistic, 
  P = res$p.value    
))


#data reduction
result$`CpG No.` <- 1:nrow(result)
result <- result[, c(ncol(result), 1:(ncol(result)-1))]

#save data
write.csv(result,file = "comparison_result.csv",row.names = F)
