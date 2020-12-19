Normal<- function (Data, Prior, Mcmc) 
{
  pandterm = function(message) {
    #print(Data)
    #print(Prior)
    #print(Mcmc)
    stop(message, call. = FALSE)
  }
  if (missing(Data)) {
    pandterm("Requires Data argument -- list of y and X")
  }
  if (is.null(Data$X)) {
    pandterm("Requires Data element X")
  }
  X = Data$X
  if (is.null(Data$y)) {
    pandterm("Requires Data element y")
  }
  y = Data$y
  nvar = ncol(X)
  nobs = length(y)
  if (nobs != nrow(X)) {
    pandterm("length(y) is not equal to nrow(X)")
  }
  if (missing(Prior)) {
    betabar = c(rep(0, nvar))
    A = 0.01 * diag(nvar)
    a = 0.001
    b = 0.001
  }
  else {
    if (is.null(Prior$betabar)) {
      betabar = c(rep(0, nvar))
    }
    else {
      betabar = Prior$betabar
    }
    if (is.null(Prior$A)) {
      A = 0.01 * diag(nvar)
    }
    else {
      A = Prior$A
    }
    if (is.null(Prior$a)) {
      a = 0.001
    }
    else {
      a = Prior$a
    }
    if (is.null(Prior$b)) {
      b = 0.001
    }
    else {
      b = Prior$b
    }
  }
  print(A)
  if (ncol(A) != nrow(A) || ncol(A) != nvar || nrow(A) != nvar) {
    pandterm(paste("bad dimensions for A", dim(A),"ojo",nvar))
  }
  if (length(betabar) != nvar) {
    pandterm(paste("betabar wrong length, length= ", length(betabar)))
  }
  if (missing(Mcmc)) {
    pandterm("requires Mcmc argument")
  }
  else {
    if (is.null(Mcmc$R)) {
      pandterm("requires Mcmc element R")
    }
    else {
      R = Mcmc$R
    }
    if (is.null(Mcmc$keep)) {
      keep = 1
    }
    else {
      keep = Mcmc$keep
    }
    if (is.null(Mcmc$burnin)) {
      burnin = 0
    }
    else {
      burnin = Mcmc$burnin
    }
    if (is.null(Mcmc$sigmasq)) {
      sigmasq = var(y)
    }
    else {
      sigmasq = Mcmc$sigmasq
    }
  }
  #   cat(" ", fill = TRUE)
  #   cat("Starting Gibbs Sampler for Univariate Regression Model", 
  #       fill = TRUE)
  #   cat("  with ", nobs, " observations", fill = TRUE)
  #   cat(" ", fill = TRUE)
  #   cat("Prior Parms: ", fill = TRUE)
  #   cat("betabar", fill = TRUE)
  #   print(betabar)
  #   cat("A", fill = TRUE)
  #   print(A)
  #   cat("nu = ", nu, " ssq= ", ssq, fill = TRUE)
  #   cat(" ", fill = TRUE)
  #   cat("MCMC parms: ", fill = TRUE)
  #   cat("R= ", R, " keep= ", keep, fill = TRUE)
  #   cat(" ", fill = TRUE)
  sigmasqdraw = double(floor((R-burnin)/keep))
  betadraw = matrix(double(floor((R-burnin) * nvar/keep)), ncol = nvar)
  XpX = crossprod(X)
  Xpy = crossprod(X, y)
  ypy = crossprod(y)
  sigmasq = as.vector(sigmasq)
  #   itime = proc.time()[3]
  #   cat("MCMC Iteration (est time to end - min) ", fill = TRUE)
  #   fsh()
  Rnew = R-burnin
  withProgress(message = 'Making calculations', value = 0, {
  for (rep in 1:Rnew) {
    B = solve(XpX + A)
    btilde = B%*%(Xpy + A%*%betabar)
    beta = mvrnorm(n = 1, btilde, sigmasq*B)
    #res = y - X %*% beta
    #s = t(res) %*% res
    sigmasq = rinvgamma(1, a + nobs, b + ypy + t(betabar)%*%A%*%betabar - t(btilde)%*%(XpX + A)%*%btilde)
    #sigmasq = (b + s)/rchisq(1, a + nobs)
    #sigmasq = as.vector(sigmasq)
    #     if (rep%%100 == 0) {
    #       ctime = proc.time()[3]
    #       timetoend = ((ctime - itime)/rep) * (R - rep)
    #       cat(" ", rep, " (", round(timetoend/60, 1), ")", 
    #           fill = TRUE)
    #       fsh()
    #     }
    
        if (rep%%keep == 0) {
      mkeep = rep/keep
      betadraw[mkeep, ] = beta
      sigmasqdraw[mkeep] = sigmasq
        }
    incProgress(1/rep, detail = paste('Doing iteration', rep))
  }
  })
  #   ctime = proc.time()[3]
  #   cat("  Total Time Elapsed: ", round((ctime - itime)/60, 2), 
  #       "\n")
  attributes(betadraw)$class = c("mcmc")
  attributes(betadraw)$mcpar = c(1, R, keep, burnin)
  attributes(sigmasqdraw)$class = c("mcmc")
  attributes(sigmasqdraw)$mcpar = c(1, R, keep, burnin)
#   betadrawsum<- summary(betadraw)
#   betadrawgeweke<- geweke.diag(betadraw[,])
#   betadrawraftery<- raftery.diag(betadraw[,],q=0.5,r=0.025,s = 0.95)
#   betadrawheidel<- heidel.diag(betadraw[,])
# plot(betadraw)
# autocorr.plot(betadraw)
#   sigmasqdrawsum<- summary(sigmasqdraw)
#   sigmasqdrawgeweke<- geweke.diag(sigmasqdraw[])
#   sigmasqdrawraftery<- raftery.diag(sigmasqdraw[],q=0.5,r=0.025,s = 0.95)
#   sigmasqdrawheidel<- heidel.diag(sigmasqdraw[])
#   return(list(betadrawsum = betadrawsum, betadrawgeweke = betadrawgeweke, betadrawraftery = betadrawraftery, betadrawheidel = betadrawheidel, sigmasqdrawsum = sigmasqdrawsum, sigmasqdrawgeweke = sigmasqdrawgeweke, sigmasqdrawraftery = sigmasqdrawraftery, sigmasqdrawheidel = sigmasqdrawheidel))
  return(list(betadraw = betadraw, sigmasqdraw = sigmasqdraw))
}

# NormalDownLoad<- function(betadraw,sigmasqdraw){
#   DownLoad<- cbind(betadraw,sigmasqdraw)
#   return(list(DownLoad))
# }


# NormalSumDiag<- function(betadraw,sigmasqdraw){
#   SummaryLocation<- summary(betadraw)
# #  names(batadrawsum)<- "Summary"
#   GewekeTestLocation<- geweke.diag(betadraw)
#   RafteryTestLocation<- raftery.diag(betadraw,q=0.5,r=0.025,s = 0.95)
#   HeidelTestLocation<- heidel.diag(betadraw)
#   SummaryScale<- summary(sigmasqdraw)
#   GewekeTestScale<- geweke.diag(sigmasqdraw)
#   RafteryTestScale<- raftery.diag(sigmasqdraw,q=0.5,r=0.025,s = 0.95)
#   HeidelTestScale<- heidel.diag(sigmasqdraw)
#   return(list(SummaryLocation = SummaryLocation, GewekeTestLocation = GewekeTestLocation, RafteryTestLocation = RafteryTestLocation, HeidelTestLocation = HeidelTestLocation, SummaryScale = SummaryScale, GewekeTestScale = GewekeTestScale, RafteryTestScale = RafteryTestScale, HeidelTestScale = HeidelTestScale))
# }



# NormalPlot<- function(betadraw){
#   hist(betadraw, breaks=20,freq=FALSE, xlab="Location", main="", col="lightgreen")
#   lines(density(betadraw,na.rm = TRUE), col="red", lwd=2)
#   abline(h = NULL, v = c(quantile(betadraw,c(0.025, 0.975))), col = "purple", lwd=2)
#   text(quantile(betadraw,c(0.025)),y=1, "Quantile 2.5%", col = "black", adj = c(0,-0.5), cex=0.75)
#   text(quantile(betadraw,c(0.975)),y=1, "Quantile 97.5%", col = "black", adj = c(0,-0.5), cex=0.75)
#   abline(h = NULL, v = c(quantile(betadraw,c(0.5))), col = "red", lwd=2)
#   abline(h = NULL, v = mean(betadraw), col = "blue", lwd=2)
#   legend("topleft",inset=.05,cex = 0.75,c("Median","Mean"),horiz=TRUE,lty=c(1,1),lwd=c(1,1),col=c("red","black"),bg="grey96")
# }


#nc<-ncol(res$betadraw[,])
#par(mfrow=c(nc,1))
#par(mar=c(1,1,1,1))
#par(mar= c(5.1, 4.1, 4.1, 2.1))
#par(mfrow=c(nc,1))
#  for (i in 1:nc ) {
#   pdf(paste("Location",paste(i,".pdf", sep = "", collapse = NULL)))
#    readline(prompt="Press [enter] to continue")
#    NormalPlot(as.vector(res$betadraw[,i]))
#   dev.off()
#   setEPS()
#   postscript(paste("Location",paste(i,".eps", sep = "", collapse = NULL)))
#   NormalPlot(as.vector(res$betadraw[,i]))
#   dev.off()
# }


# path<-getwd()
#setwd(paste(path,"/Insumo",sep=""))
# setwd("..")
# dir.create(file.path(path,"Posteriors"),showWarnings = FALSE)
# setwd(file.path(path,"Posteriors"))


# NormalPlot<- function(betadraw,sigmasqdraw){
#   n<- ncol(betadraw)+1
#   par(mfrow=c(n, 2))
#   plot(betadraw)
# #  autocorr.plot(betadraw)
#   plot(sigmasqdraw)
# #  autocorr.plot(sigmasqdraw)
# }

# setwd("C:/andres_ramirez/Andres_EAFIT/Bayesian Econometrics/UserInterface/BEsmartV18/Data")
# Data1<- read.table(file="Datalinear.csv",header=TRUE,sep=",")
# Data<-list(y=as.vector(Data1[,1]),X=as.matrix(Data1[,-1]))
# Bmean<-c(0,0,0,0)
# Bvar<- 0.01*diag(4)
# a<- 3
# b<- var(Data1[,1])
# Prior<- list(Bmean,Bvar,a,b)
# MCMC<- list(R=10000,keep=20,burnin=100)
# 
# res<-Normal(Data,Prior,MCMC)
# NormalSumDiag(res$betadraw[,],res$sigmasq[])





