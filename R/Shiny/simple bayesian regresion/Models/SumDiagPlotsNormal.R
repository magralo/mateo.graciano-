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

# SumDiag<- function(betadraw){
#   SummaryLocation<- summary(betadraw)
#   GewekeTestLocation<- geweke.diag(betadraw)
#   RafteryTestLocation<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
#   HeidelTestLocation<- heidel.diag(betadraw)
#   return(list(SummaryLocation = SummaryLocation, GewekeTestLocation = GewekeTestLocation, RafteryTestLocation = RafteryTestLocation, HeidelTestLocation = HeidelTestLocation))
# }

# Plot<- function(betadraw){
#   hist(betadraw, breaks=20,freq=FALSE, xlab="Parameter", main="", col="lightgreen")
#   lines(density(betadraw,na.rm = TRUE), col="red", lwd=2)
#   abline(h = NULL, v = c(quantile(betadraw,c(0.025, 0.975))), col = "purple", lwd=2)
#   text(quantile(betadraw,c(0.025)),y=1, "Quantile 2.5%", col = "black", adj = c(0,-0.5), cex=0.75)
#   text(quantile(betadraw,c(0.975)),y=1, "Quantile 97.5%", col = "black", adj = c(0,-0.5), cex=0.75)
#   abline(h = NULL, v = c(quantile(betadraw,c(0.5))), col = "red", lwd=2)
#   abline(h = NULL, v = mean(betadraw), col = "blue", lwd=2)
#   legend("topleft",inset=.05,cex = 0.75,c("Median","Mean"),horiz=TRUE,lty=c(1,1),lwd=c(1,1),col=c("red","black"),bg="grey96")
# }
# 
# Plot.trace<- function(betadraw){traceplot(mcmc(betadraw), main = "Trace Plot", xlab="Iteration", ylab= "Parameter", col= "blue")}
# Plot.corr<- function(betadraw){autocorr.plot(mcmc(betadraw), main = "Autocorrelation", col="blue")}

