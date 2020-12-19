# SumDiag<- function(betadraw,sigmasqdraw){
#   SummaryLocation<- summary(betadraw)
#   #  names(batadrawsum)<- "Summary"
#   GewekeTestLocation<- geweke.diag(betadraw)
#   RafteryTestLocation<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
#   HeidelTestLocation<- heidel.diag(betadraw)
#   SummaryScale<- summary(sigmasqdraw)
#   GewekeTestScale<- geweke.diag(sigmasqdraw)
#   RafteryTestScale<- raftery.diag(sigmasqdraw,q=0.5,r=0.025,s = 0.95)
#   HeidelTestScale<- heidel.diag(sigmasqdraw)
#   return(list(SummaryLocation = SummaryLocation, GewekeTestLocation = GewekeTestLocation, RafteryTestLocation = RafteryTestLocation, HeidelTestLocation = HeidelTestLocation, SummaryScale = SummaryScale, GewekeTestScale = GewekeTestScale, RafteryTestScale = RafteryTestScale, HeidelTestScale = HeidelTestScale))
# }

model.formula<- function(formula, data=list(), ...)
{
  mf <- model.frame(formula=formula, data=data)
  X <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  return(list(X=X,y=y))
}

model.formula1<- function(formula, data=list(), ...)
{
  mf <- model.frame(formula=formula, data=data)
  Z <- model.matrix(attr(mf, "terms"), data=mf)
  x <- model.response(mf)
  return(list(Z=Z,x=x))
}

SumDiagNormal<- function(betadraw,sigmasqdraw){
  SummaryLocation<- summary(betadraw)
  GewekeTestLocation<- geweke.diag(betadraw)
  RafteryTestLocation<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLocation<- heidel.diag(betadraw)
  SummaryScale<- summary(sigmasqdraw)
  GewekeTestScale<- geweke.diag(sigmasqdraw)
  RafteryTestScale<- raftery.diag(sigmasqdraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestScale<- heidel.diag(sigmasqdraw)
  return(list(SummaryLocation = SummaryLocation, GewekeTestLocation = GewekeTestLocation, RafteryTestLocation = RafteryTestLocation, HeidelTestLocation = HeidelTestLocation, SummaryScale = SummaryScale, GewekeTestScale = GewekeTestScale, RafteryTestScale = RafteryTestScale, HeidelTestScale = HeidelTestScale))
}

SumDiagLogit<- function(betadraw){
  Summary<- summary(betadraw)
  GewekeTest<- geweke.diag(betadraw)
  RafteryTest<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTest<- heidel.diag(betadraw)
  return(list(Summary = Summary, Geweke.Test = GewekeTest, Raftery.Test = RafteryTest, Heidel.Test = HeidelTest))
}

SumDiagProbit<- function(betadraw){
  Summary<- summary(betadraw)
  GewekeTest<- geweke.diag(betadraw)
  RafteryTest<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTest<- heidel.diag(betadraw)
  return(list(Summary = Summary, Geweke.Test = GewekeTest, Raftery.Test = RafteryTest, Heidel.Test = HeidelTest))
}

SumDiagMultProbit<- function(betadraw,sigmadraw){
  SummaryDependentVar<- summary(betadraw)
  GewekeTestDependentVar<- geweke.diag(betadraw)
  RafteryTestDependentVar<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestDependentVar<- heidel.diag(betadraw)
  SummaryCovMatrix<- summary(sigmadraw)
  GewekeTestCovMatrix<- geweke.diag(sigmadraw)
  RafteryTestCovMatrix<- raftery.diag(sigmadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestCovMatrix<- heidel.diag(sigmadraw)
  return(list(SummaryDependentVar = SummaryDependentVar, GewekeTestDependentVar = GewekeTestDependentVar, RafteryTestDependentVar = RafteryTestDependentVar, HeidelTestDependentVar = HeidelTestDependentVar, SummaryCovMatrix = SummaryCovMatrix, GewekeTestCovMatrix = GewekeTestCovMatrix, RafteryTestCovMatrix = RafteryTestCovMatrix, HeidelTestCovMatrix = HeidelTestCovMatrix))
}

SumDiagMultLogit<- function(betadraw){
  Summary<- summary(betadraw)
  GewekeTest<- geweke.diag(betadraw)
  RafteryTest<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTest<- heidel.diag(betadraw)
  return(list(Summary = Summary, Geweke.Test = GewekeTest, Raftery.Test = RafteryTest, Heidel.Test = HeidelTest))
}


SumDiagOprobit<- function(betadraw,cutdraw){
  SummaryLocation<- summary(betadraw)
  GewekeTestLocation<- geweke.diag(betadraw)
  RafteryTestLocation<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLocation<- heidel.diag(betadraw)
  SummaryCut<- summary(cutdraw)
  GewekeTestCut<- geweke.diag(cutdraw)
  RafteryTestCut<- raftery.diag(cutdraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestCut<- heidel.diag(cutdraw)
  return(list(SummaryLocation = SummaryLocation, GewekeTestLocation = GewekeTestLocation, RafteryTestLocation = RafteryTestLocation, HeidelTestLocation = HeidelTestLocation, SummaryCut = SummaryCut, GewekeTestCut = GewekeTestCut, RafteryTestCut = RafteryTestCut, HeidelTestCut = HeidelTestCut))
}

SumDiagNegBin<- function(betadraw,alphadraw){
  SummaryLocation<- summary(betadraw)
  GewekeTestLocation<- geweke.diag(betadraw)
  RafteryTestLocation<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLocation<- heidel.diag(betadraw)
  SummaryDispersion<- summary(alphadraw)
  GewekeTestDispersion<- geweke.diag(alphadraw)
  RafteryTestDispersion<- raftery.diag(alphadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestDispersion<- heidel.diag(alphadraw)
  return(list(SummaryLocation = SummaryLocation, GewekeTestLocation = GewekeTestLocation, RafteryTestLocation = RafteryTestLocation, HeidelTestLocation = HeidelTestLocation, SummaryDispersion = SummaryDispersion, GewekeTestDispersion = GewekeTestDispersion, RafteryTestDispersion = RafteryTestDispersion, HeidelTestDispersion = HeidelTestDispersion))
}

SumDiagTobit<- function(betadraw){
  Summary<- summary(betadraw)
  GewekeTest<- geweke.diag(betadraw)
  RafteryTest<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTest<- heidel.diag(betadraw)
  return(list(Summary = Summary, Geweke.Test = GewekeTest, Raftery.Test = RafteryTest, Heidel.Test = HeidelTest))
}

SumDiagQuantile<- function(betadraw){
  Summary<- summary(betadraw)
  GewekeTest<- geweke.diag(betadraw)
  RafteryTest<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTest<- heidel.diag(betadraw)
  return(list(Summary = Summary, Geweke.Test = GewekeTest, Raftery.Test = RafteryTest, Heidel.Test = HeidelTest))
}

SumDiagHier<- function(betadraw){
  Summary<- summary(betadraw)
  GewekeTest<- geweke.diag(betadraw)
  RafteryTest<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTest<- heidel.diag(betadraw)
  return(list(Summary = Summary, Geweke.Test = GewekeTest, Raftery.Test = RafteryTest, Heidel.Test = HeidelTest))
}

SumDiagInstVar<- function(betadraw,deltadraw,gammadraw,Sigmadraw){
  SummaryEndogenousVar<- summary(betadraw)
  GewekeTestEndogenousVar<- geweke.diag(betadraw)
  RafteryTestEndogenousVar<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestEndogenousVar<- heidel.diag(betadraw)
  SummaryInstruments<- summary(deltadraw)
  GewekeTestInstruments<- geweke.diag(deltadraw)
  RafteryTestInstruments<- raftery.diag(deltadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestInstruments<- heidel.diag(deltadraw)
  SummaryExogenousVars<- summary(gammadraw)
  GewekeTestExogenousVars<- geweke.diag(gammadraw)
  RafteryTestExogenousVars<- raftery.diag(gammadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestExogenousVars<- heidel.diag(gammadraw)
  SummaryCovMatrix<- summary(Sigmadraw)
  GewekeTestCovMatrix<- geweke.diag(Sigmadraw)
  RafteryTestCovMatrix<- raftery.diag(Sigmadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestCovMatrix<- heidel.diag(Sigmadraw)
  return(list(SummaryEndogenousVar = SummaryEndogenousVar, GewekeTestEndogenousVar = GewekeTestEndogenousVar, RafteryTestEndogenousVar = RafteryTestEndogenousVar, HeidelTestEndogenousVar = HeidelTestEndogenousVar, SummaryInstruments = SummaryInstruments, GewekeTestInstruments = GewekeTestInstruments, RafteryTestInstruments = RafteryTestInstruments, HeidelTestInstruments = HeidelTestInstruments,
              SummaryExogenousVars = SummaryExogenousVars, GewekeTestExogenousVars = GewekeTestExogenousVars, RafteryTestExogenousVars = RafteryTestExogenousVars, HeidelTestExogenousVars = HeidelTestExogenousVars, SummaryCovMatrix = SummaryCovMatrix, GewekeTestCovMatrix = GewekeTestCovMatrix, RafteryTestCovMatrix = RafteryTestCovMatrix, HeidelTestCovMatrix = HeidelTestCovMatrix))
}

SumDiagMultiReg<- function(betadraw,sigmasqdraw){
  SummaryLocation<- summary(betadraw)
  GewekeTestLocation<- geweke.diag(betadraw)
  RafteryTestLocation<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLocation<- heidel.diag(betadraw)
  SummaryScale<- summary(sigmasqdraw)
  GewekeTestScale<- geweke.diag(sigmasqdraw)
  RafteryTestScale<- raftery.diag(sigmasqdraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestScale<- heidel.diag(sigmasqdraw)
  return(list(SummaryLocation = SummaryLocation, GewekeTestLocation = GewekeTestLocation, RafteryTestLocation = RafteryTestLocation, HeidelTestLocation = HeidelTestLocation, SummaryScale = SummaryScale, GewekeTestScale = GewekeTestScale, RafteryTestScale = RafteryTestScale, HeidelTestScale = HeidelTestScale))
}

SumDiagSUR<- function(betadraw,sigmasqdraw){
  SummaryLocation<-summary(betadraw)
  GewekeTestLocation<- geweke.diag(mcmc(betadraw))
  RafteryTestLocation<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLocation<- heidel.diag(betadraw)
  SummaryScale<- summary(sigmasqdraw)
  GewekeTestScale<- geweke.diag(mcmc(sigmasqdraw))
  RafteryTestScale<- raftery.diag(sigmasqdraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestScale<- heidel.diag(sigmasqdraw)
  return(list(SummaryLocation = SummaryLocation, GewekeTestLocation = GewekeTestLocation, RafteryTestLocation = RafteryTestLocation, HeidelTestLocation = HeidelTestLocation, SummaryScale = SummaryScale, GewekeTestScale = GewekeTestScale, RafteryTestScale = RafteryTestScale, HeidelTestScale = HeidelTestScale))
}


SumDiagBVProbit<- function(betadraw,rhodraw){
  SummaryLocation<-summary(betadraw)
  GewekeTestLocation<- geweke.diag(mcmc(betadraw))
  RafteryTestLocation<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestLocation<- heidel.diag(betadraw)
  SummaryCorrelation<- summary(rhodraw)
  GewekeTestCorrelation<- geweke.diag(mcmc(rhodraw))
  RafteryTestCorrelation<- raftery.diag(rhodraw,q=0.5,r=0.025,s = 0.95)
  HeidelTestCorrelation<- heidel.diag(rhodraw)
  return(list(SummaryLocation = SummaryLocation, GewekeTestLocation = GewekeTestLocation, RafteryTestLocation = RafteryTestLocation, HeidelTestLocation = HeidelTestLocation, SummaryCorrelation = SummaryCorrelation, GewekeTestCorrelation = GewekeTestCorrelation, RafteryTestCorrelation = RafteryTestCorrelation, HeidelTestCorrelation = HeidelTestCorrelation))
}

SumDiagBayBoots<- function(betadraw){
  Summary<- summary(betadraw)
  GewekeTest<- geweke.diag(betadraw)
  RafteryTest<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
  HeidelTest<- heidel.diag(betadraw)
  return(list(Summary = Summary, Geweke.Test = GewekeTest, Raftery.Test = RafteryTest, Heidel.Test = HeidelTest))
}


Plot<- function(betadraw){
  hist(betadraw, breaks=20,freq=FALSE, xlab="Parameter", main="", col="lightgreen")
  lines(density(betadraw,na.rm = TRUE), col="red", lwd=2)
  abline(h = NULL, v = c(quantile(betadraw,c(0.025, 0.975))), col = "purple", lwd=2)
  text(quantile(betadraw,c(0.025)),y=1, "Quantile 2.5%", col = "black", adj = c(0,-0.5), cex=0.75)
  text(quantile(betadraw,c(0.975)),y=1, "Quantile 97.5%", col = "black", adj = c(0,-0.5), cex=0.75)
  abline(h = NULL, v = c(quantile(betadraw,c(0.5))), col = "red", lwd=3)
  abline(h = NULL, v = mean(betadraw), col = "blue", lwd=2)
  legend("topleft",inset=.05,cex = 0.75,c("Median","Mean"),horiz=TRUE,lty=c(1,1),lwd=c(1,1),col=c("red","blue"),bg="grey96")
}

Plot.trace<- function(betadraw){traceplot(mcmc(betadraw), main = "Trace Plot", xlab="Iteration", ylab= "Parameter", col= "blue")}
Plot.corr<- function(betadraw){autocorr.plot(mcmc(betadraw), main = "Autocorrelation", col="blue")}
