

################################################################################
# Functions to conduct the three ICE FALCON models
################################################################################


################################################################################
# Function to convert the data from one twin to other twin
################################################################################

convertfn<-function(data)
{
   x<-data
   n<-length(x);

   twin1seq<-seq(1,n,2);                        # generate sequence id for twin 1
   twin2seq<-seq(2,n,2);                        # generate sequence id for twin 2

   xtwin1<-x;
   xtwin2<-x;

# Assign twin 1 the value of 1  and reverse the order of twin 2

   xtwin1[twin1seq]<-1;
   xtwin11<-xtwin1[-1];
   twin1<-c(xtwin11,1);

# Assign twin 2 the value of 1 and reverse the order of twin 1

   xtwin2[twin2seq]<-1;
   xtwin21<-c(1,xtwin2);
   twin2<-xtwin21[-length(xtwin21)];

# Compute co-twin by element-by-element multiplication of vector

   cotwinx<-twin1*twin2;
   return(cotwinx)
}



################################################################################
# Function to compute estimate regression coefficients for three models in
# the ICE FALCON using exchangeable working correlation structure in GEE
################################################################################

gee.coef<-function(data){

             gee.fit1<-geeglm(mod1.formula,data=data,id=familyid,corstr="exchangeable")           # Model I
             gee.fit2<-geeglm(mod2.formula,data=data,id=familyid,corstr="exchangeable")           # Model II
             gee.fit3<-geeglm(mod3.formula,data=data,id=familyid,corstr="exchangeable")           # Model III

             coef1 <- gee.fit1$coef[2]
             coef2 <- gee.fit2$coef[2]
             coef3 <- gee.fit3$coef[2]
             coef4 <- gee.fit3$coef[3]

             list(betaself=coef1,betacotwin=coef2,betaself.ad=coef3,betacotwin.ad=coef4)
}



#####################################################################################
# Function to generate bootstrap data and compute bootstrap estimates for the three
# ICE FALCON model using GEE
#####################################################################################


gee.boot<-function(data,m){                                            # m = number of boostraps

       dat<-data[order(data$twinid,data$familyid),]                    # sort the data according to twinid

       coef1<-c(0)
       coef2<-c(0)
       coef3<-c(0)
       coef4<-c(0)

      for(i in 1:m){
       n<-nrow(dat)/2                             # n = the number of twin pairs
       samp<-sample(1:n,replace=T)
       rep<-rep(n,n)+samp
       com<-c(samp,rep)
       inputdat<-dat[com,]
       id2<-rep(1:n,2)
       inputdat2<-cbind(id2,inputdat)
       inputdat3<-inputdat2[order(inputdat2$id2),]
       # print(id2)
       # print(inputdat3)
       mod1<-geeglm(mod1.formula,data=inputdat3,id=id2,corstr="exchangeable")
       mod2<-geeglm(mod2.formula,data=inputdat3,id=id2,corstr="exchangeable")
       mod3<-geeglm(mod3.formula,data=inputdat3,id=id2,corstr="exchangeable")
       coef1[i]<-coefficients(mod1)[2]
       coef2[i]<-coefficients(mod2)[2]
       coef3[i]<-coefficients(mod3)[2]
       coef4[i]<-coefficients(mod3)[3]
      }
      list(betaself=coef1,betacotwin=coef2,betaself.ad=coef3,betacotwin.ad=coef4)
}




#################################################################################
# Function to compute change in coefficient and its p-value
#################################################################################

coef.change<-function(coef.model,coef.boot){

       coef.diff.self<-coef.model$betaself.ad-coef.model$betaself
       coef.diff.cot <-coef.model$betacotwin.ad-coef.model$betacotwin

       diff.self<-coef.boot$betaself.ad-coef.boot$betaself
       diff.cotwin<-coef.boot$betacotwin.ad-coef.boot$betacotwin

       se.self<- sd(diff.self)
       se.cotwin<-sd(diff.cotwin)

       z.self<-coef.diff.self/se.self
       z.cotwin<-coef.diff.cot/se.cotwin
       p.self<-2*pnorm(abs(z.self),lower.tail=F)
       p.cotwin<-2*pnorm(abs(z.cotwin),lower.tail=F)

       list(beta_self_change = coef.diff.self, SE_self_change = se.self, pvalue_self_change = p.self, 
            beta_cot_change =  coef.diff.cot, SE_cotwin_change = se.cotwin, pvalue_cotwin_change = p.cotwin)
}


