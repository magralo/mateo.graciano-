MultinomialLogit<-function (Data, Prior, Mcmc) 
{
  library(MCMCpack)
  pandterm = function(message) {
    stop(message, call. = FALSE)
  }
  if (missing(Data)) {
    pandterm("Requires Data argument -- list of p, y, X")
  }
  if (is.null(Data$X)) {
    pandterm("Requires Data element X")
  }
  X = Data$X
  if (is.null(Data$y)) {
    pandterm("Requires Data element y")
  }
  y = Data$y
  if (is.null(Data$p)) {
    pandterm("Requires Data element p")
  }
  p = Data$p
  nvar = ncol(X)
  nobs = length(y)
  if (length(y) != (nrow(X)/p)) {
    pandterm("length(y) ne nrow(X)/p")
  }
  if (sum(y %in% (1:p)) < nobs) {
    pandterm("invalid values in y vector -- must be integers in 1:p")
  }
#   cat(" table of y values", fill = TRUE)
#   print(table(y))
  if (missing(Prior)) {
    betabar = c(rep(0, nvar))
    A = 0.01 * diag(nvar)
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
  }
  if (ncol(A) != nrow(A) || ncol(A) != nvar || nrow(A) != nvar) {
    pandterm(paste("bad dimensions for A", dim(A)))
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
    if (is.null(Mcmc$nu)) {
      nu = 6
    }
    else {
      nu = Mcmc$nu
    }
  }
#   cat(" ", fill = TRUE)
#   cat("Starting Independence Metropolis Sampler for Multinomial Logit Model", 
#       fill = TRUE)
#   cat("  ", length(y), " obs with ", p, " alternatives", fill = TRUE)
#   cat(" ", fill = TRUE)
#   cat("Table of y Values", fill = TRUE)
#   print(table(y))
#   cat("Prior Parms: ", fill = TRUE)
#   cat("betabar", fill = TRUE)
#   print(betabar)
#   cat("A", fill = TRUE)
#   print(A)
#   cat(" ", fill = TRUE)
#   cat("MCMC parms: ", fill = TRUE)
#   cat("R= ", R, " keep= ", keep, " nu (df for st candidates) = ", 
#       nu, fill = TRUE)
#   cat(" ", fill = TRUE)
  Rnew = R-burnin
  betadraw = matrix(double(floor(Rnew/keep) * nvar), ncol = nvar)
  loglike = double(floor(Rnew/keep))
  beta = c(rep(0, nvar))
  mle = optim(beta, llmnl, X = X, y = y, method = "BFGS", hessian = TRUE, 
              control = list(fnscale = -1))
  beta = mle$par
  betastar = mle$par
  mhess = mnlHess(beta, y, X)
  candcov = chol2inv(chol(mhess))
  root = chol(candcov)
  rooti = backsolve(root, diag(nvar))
  priorcov = chol2inv(chol(A))
  rootp = chol(priorcov)
  rootpi = backsolve(rootp, diag(nvar))
#   itime = proc.time()[3]
#   cat("MCMC Iteration (est time to end - min) ", fill = TRUE)
#   fsh()
  oldloglike = llmnl(beta, y, X)
  oldlpost = oldloglike + lndMvn(beta, betabar, rootpi)
  oldlimp = lndMvst(beta, nu, betastar, rooti)
  naccept = 0
  
  ## ITER
  withProgress(message = 'Making calculations', value = 0, {
  for (rep in 1:Rnew) {
    betac = rmvst(nu, betastar, root)
    cloglike = llmnl(betac, y, X)
    clpost = cloglike + lndMvn(betac, betabar, rootpi)
    climp = lndMvst(betac, nu, betastar, rooti)
    ldiff = clpost + oldlimp - oldlpost - climp
    alpha = min(1, exp(ldiff))
    if (alpha < 1) {
      unif = runif(1)
    }
    else {
      unif = 0
    }
    if (unif <= alpha) {
      beta = betac
      oldloglike = cloglike
      oldlpost = clpost
      oldlimp = climp
      naccept = naccept + 1
    }
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
      loglike[mkeep] = oldloglike
    }
    incProgress(1/rep, detail = paste('Doing iteration', rep))
  }
  })
#   ctime = proc.time()[3]
#   cat("  Total Time Elapsed: ", round((ctime - itime)/60, 2), 
#       "\n")
#   attributes(betadraw)$class = c("bayesm.mat", "mcmc")
  attributes(betadraw)$class = c("mcmc")
  attributes(betadraw)$mcpar = c(1, R, keep, burnin)
#   return(list(betadraw = betadraw, loglike = loglike, acceptr = naccept/R))
return(list(betadraw = betadraw))
}


# Data<- read.table(file="SimMultLogitmodel.csv",header=TRUE,sep=",")
# Xa<-as.matrix(Data[,2:10])
# Xd<-as.matrix(Data[,11])
# XMPP<- createX(3, na=3, nd=1, Xa=Xa, Xd=Xd, INT = TRUE, DIFF = FALSE, base = 3)
# #simout=simmnp(X,p,500,beta,Sigma)
# Data1=list(y=as.vector(Data[,1]),X=XMPP,p=3)
# Mcmc1=list(R=1000,keep=1)
# out=MultLogit(Data=Data1,Mcmc=Mcmc1)
# summary(out$betadraw)
