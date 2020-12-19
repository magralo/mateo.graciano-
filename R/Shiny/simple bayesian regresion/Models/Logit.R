remove(list = ls())
setwd("C:/andres_ramirez/Andres_EAFIT/Bayesian Econometrics/UserInterface/BEsmarterV1/Data")
Data1<- read.table(file="SimLogitmodel.csv",header=TRUE,sep=",")
#attach(Data1)
Data<-list(y=as.vector(Data1[,1]),X=as.matrix(Data1[,-1]))
form<-y~x1+x2
Bmean<-c(0,0,0)
Bvar<- 0.01*diag(3)
MCMC<-list(mcmc=11000,burnin=1000,thin=10)
list1<-list(form, data=Data, burnin = 1000, mcmc = 10000,thin=1, tune=1.1, verbose = 0, seed = NA, beta.start = NA, b0 = Bmean, B0 = Bvar)
#list<-list(form,data=Data,b0=Bmean,B0=Bvar,mcmc=11000,burnin=1000,thin=10)
#res<-MCMClogit(form=list1[[1]], data=list1[[2]], burnin = list1[[3]], mcmc = list1[[4]],thin=list1[[5]], tune=1.1, verbose = 0, seed = NA, beta.start = NA, b0 = list1[[10]], B0 = list1[[11]])
#summary(res)
res1<-do.call(MCMClogit,list1)
summary(res1)



