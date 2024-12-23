rm(list = ls()) 

#load data
phen_file <- "sample_information1.RDS"
meth_file <- "top_CpG_level.RDS"
result_file <- "eGFR.RDS"
plot_file <- "plot.pdf"

phen <- readRDS(phen_file)
meth <- readRDS(meth_file)

if(!all(colnames(meth) == phen$no)) stop("mismatch id")
if(!all(table(phen$family_ID) == 2)) stop("not all are twins")

phen<-phen[order(phen$family_ID,phen$no),]

phen$family_ID <- as.character(phen$family_ID)

phen$no <- as.character(phen$no)

meth <- meth[, phen$no]

########### cause inference ##############
phen$phen <- phen$eGFR

dat <- phen

dat$outcome <- dat$phen

## x2y methylation -> phenotype
result <- list()
for(i in 1:nrow(meth)) {
  
  dat$predictor <- as.numeric(meth[i, ])
  unitname <- rownames(meth)[i]
  
  dat$familyid <- dat$family_ID
  dat$twinid <- dat$no
  
  library(geepack)
  
  source("ICEFALCONfunction.r")  
  
  dat$y<-dat$outcome
  dat$x<-dat$predictor
  dat$ycot <- 0
  dat$xcot <- 0
  dat$ycot[dat$intrapair == 1] <- dat$y[dat$intrapair == 2]
  dat$ycot[dat$intrapair == 2] <- dat$y[dat$intrapair == 1]
  dat$xcot[dat$intrapair == 1] <- dat$x[dat$intrapair == 2]
  dat$xcot[dat$intrapair == 2] <- dat$x[dat$intrapair == 1]
  
  mod1.formula<-as.formula( "y~x")
  mod2.formula<-as.formula( "y~xcot")
  mod3.formula<-as.formula( "y~x+xcot")
  
  gee.fit1<-geeglm(mod1.formula,data=dat,id=familyid,corstr="exchangeable")
  gee.fit2<-geeglm(mod2.formula,data=dat,id=familyid,corstr="exchangeable")
  gee.fit3<-geeglm(mod3.formula,data=dat,id=familyid,corstr="exchangeable")  
  
  coef.model.gee<-gee.coef(dat)
  
  # Compute bootstrap estimates for the regression coefficients
  set.seed(1)
  m<-1000                                              # m = 1000 bootstraps
  coef.boot.gee<-gee.boot(dat,m)
  
  # Compute parameters change and their p-value
  
  out.gee<-coef.change(coef.model.gee,coef.boot.gee)
  # print(out.gee)
  rst <- data.frame(unit = unitname, out.gee)
  
  result[[i]] <- rst
}

#save result
require(data.table)
result <- as.data.frame(rbindlist(result))
colnames(result)[2:ncol(result)] <- paste0("x2y_", colnames(result)[2:ncol(result)])
result_x2y <- result


## x2y phenotype -> methylation
result <- list()
for(i in 1:nrow(meth)) {
  
  dat$predictor <- as.numeric(meth[i, ])
  unitname <- rownames(meth)[i]
  
  dat$familyid <- dat$family_ID
  dat$twinid <- dat$no
  
  library(geepack)
  
  source("ICEFALCONfunction.r")
  
  dat$y<-dat$predictor
  dat$x<-dat$outcome
  dat$ycot <- 0
  dat$xcot <- 0
  dat$ycot[dat$intrapair == 1] <- dat$y[dat$intrapair == 2]
  dat$ycot[dat$intrapair == 2] <- dat$y[dat$intrapair == 1]
  dat$xcot[dat$intrapair == 1] <- dat$x[dat$intrapair == 2]
  dat$xcot[dat$intrapair == 2] <- dat$x[dat$intrapair == 1]
  
  mod1.formula<-as.formula( "y~x")
  mod2.formula<-as.formula( "y~xcot")
  mod3.formula<-as.formula( "y~x+xcot")
  
  gee.fit1<-geeglm(mod1.formula,data=dat,id=familyid,corstr="exchangeable")
  gee.fit2<-geeglm(mod2.formula,data=dat,id=familyid,corstr="exchangeable")
  gee.fit3<-geeglm(mod3.formula,data=dat,id=familyid,corstr="exchangeable")  
  
  coef.model.gee<-gee.coef(dat)
  
  # Compute bootstrap estimates for the regression coefficients
  set.seed(1)
  m<-1000                                              # m = 1000 bootstraps
  coef.boot.gee<-gee.boot(dat,m)
  
  # Compute parameters change and their p-value
  
  out.gee<-coef.change(coef.model.gee,coef.boot.gee)
  # print(out.gee)
  rst <- data.frame(unit = unitname, out.gee)
  
  result[[i]] <- rst
}

#save result
require(data.table)
result <- as.data.frame(rbindlist(result))
colnames(result)[2:ncol(result)] <- paste0("y2x_", colnames(result)[2:ncol(result)])
result_y2x <- result

#merge result
require(plyr)
result <- join(result_x2y, result_y2x, type = "full")
saveRDS(result, result_file)

#save final result
saveRDS(result, "result.RDS")
write.csv(result, "result.csv")


#change the creas for CpG sites
myvars <- row.names(result)
result$myvars <- myvars
library(plyr)
result <- rename(result, c(unit="unit_former", myvars="unit"))


#plot
require(ggplot2)
require(ggrepel)
require(cowplot)

limit <- max(abs(c(result$x2y_beta_cot_change, result$x2y_SE_self_change))) * 1.1
plota <- ggplot(data = result, mapping = aes(x = x2y_beta_cot_change, y = x2y_SE_self_change, label = unit)) +
  geom_hline(yintercept = 0, linetype = "dotdash") +
  geom_vline(xintercept = 0, linetype = "dotdash") +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(colour = "red") +
  geom_text_repel(colour = "red") +
  coord_equal(xlim=c(-limit, limit),ylim=c(-limit, limit)) +
  xlab("beta co-twin change (methylation -> eGFR)") +
  ylab("beta self change (methylation -> eGFR)") +
  theme_bw()


limit <- max(abs(c(result$y2x_beta_cot_change, result$y2x_SE_self_change))) * 1.1
plotb <- ggplot(data = result, mapping = aes(x = y2x_beta_cot_change, y = y2x_SE_self_change, label = unit)) +
  geom_hline(yintercept = 0, linetype = "dotdash") +
  geom_vline(xintercept = 0, linetype = "dotdash") +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(colour = "red") +
  geom_text_repel(colour = "red") +
  coord_equal(xlim=c(-limit, limit),ylim=c(-limit, limit)) +
  xlab("beta co-twin change (eGFR -> methylation)") +
  ylab("beta self change (eGFR -> methylation)") +
  theme_bw()

pdf(plot_file, paper = "a4r", width = 0, height = 0)
plot_grid(plota, plotb, labels = c("A", "B"))
dev.off()






