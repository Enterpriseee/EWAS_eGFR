#load data
data<-fread("data_with_MRS.csv")

#analysis
eGFR<-data$eGFR
corc<-as.factor(data$corc)
MRS<-data$MRS
age<-data$age
gender<-as.factor(data$gender)
UA<-data$UA
GLU<-data$GLU
CHOL<-data$CHOL
TG<-data$TG
HDLC<-data$HDLC

result <- glm(corc ~ MRS + age + gender + UA + GLU + TG + CHOL +HDLC, data = data, family = binomial)
summary(result)
capture.output(result,file = "result_MRS_corc.csv")
