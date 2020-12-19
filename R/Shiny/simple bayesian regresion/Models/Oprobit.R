Oprobit<- function (Data, Prior, Mcmc) 
{
  breg1 = function(root, X, y, Abetabar) {
    cov = crossprod(root, root)
    betatilde = cov %*% (crossprod(X, y) + Abetabar)
    betatilde + t(root) %*% rnorm(length(betatilde))
  }
  dstartoc = function(dstar) {
    c(-100, 0, cumsum(exp(dstar)), 100)
  }
  lldstar = function(dstar, y, mu) {
    gamma = dstartoc(dstar)
    arg = pnorm(gamma[y + 1] - mu) - pnorm(gamma[y] - mu)
    epsilon = 1e-50
    arg = ifelse(arg < epsilon, epsilon, arg)
    return(sum(log(arg)))
  }
  dstarRwMetrop = function(y, mu, olddstar, s, inc.root, dstarbar, 
                           oldll, rootdi) { ### esta s es la que debemos poner
    stay = 0
    dstarc = olddstar + s * t(inc.root) %*% (matrix(rnorm(ncut), 
                                                    ncol = 1))
    cll = lldstar(dstarc, y, mu)
    clpost = cll + lndMvn(dstarc, dstarbar, rootdi)
    ldiff = clpost - oldll - lndMvn(olddstar, dstarbar, rootdi)
    alpha = min(1, exp(ldiff))
    if (alpha < 1) {
      unif = runif(1)
    }
    else {
      unif = 0
    }
    if (unif <= alpha) {
      dstardraw = dstarc
      oldll = cll
    }
    else {
      dstardraw = olddstar
      stay = 1
    }
    return(list(dstardraw = dstardraw, oldll = oldll, stay = stay))
  }
  pandterm = function(message) {
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
#   if (is.null(Data$k)) {
#     pandterm("Requires Data element k")
#   }
#  k = Data$k
  k=max(y)
  nvar = ncol(X)
  nobs = length(y)
  ndstar = k - 2
  ncuts = k + 1
  ncut = ncuts - 3
  if (length(y) != nrow(X)) {
    pandterm("y and X not of same row dim")
  }
#   if (sum(unique(y) %in% (1:k)) < length(unique(y))) {
#     pandterm("some value of y is not vaild")
#   }
  if (missing(Prior)) {
    betabar = c(rep(0, nvar))
    A = 0.01 * diag(nvar)
    Ad = diag(ndstar)
    dstarbar = c(rep(0, ndstar))
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
    if (is.null(Prior$Ad)) {
      Ad = diag(ndstar)
    }
    else {
      Ad = Prior$Ad
    }
    if (is.null(Prior$dstarbar)) {
      dstarbar = c(rep(0, ndstar))
    }
    else {
      dstarbar = Prior$dstarbar
    }
  }
  if (ncol(A) != nrow(A) || ncol(A) != nvar || nrow(A) != nvar) {
    pandterm(paste("bad dimensions for A", dim(A)))
  }
  if (length(betabar) != nvar) {
    pandterm(paste("betabar wrong length, length= ", length(betabar)))
  }
  if (ncol(Ad) != nrow(Ad) || ncol(Ad) != ndstar || nrow(Ad) != 
        ndstar) {
    pandterm(paste("bad dimensions for Ad", dim(Ad)))
  }
  if (length(dstarbar) != ndstar) {
    pandterm(paste("dstarbar wrong length, length= ", length(dstarbar)))
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
    if (is.null(Mcmc$s)) {
      s = 2.93/sqrt(ndstar)
    }
    else {
      s = Mcmc$s
    }
  }
#   cat(" ", fill = TRUE)
#   cat("Starting Gibbs Sampler for Ordered Probit Model", fill = TRUE)
#   cat("   with ", nobs, "observations", fill = TRUE)
#   cat(" ", fill = TRUE)
#   cat("Table of y values", fill = TRUE)
#   print(table(y))
#   cat(" ", fill = TRUE)
#   cat("Prior Parms: ", fill = TRUE)
#   cat("betabar", fill = TRUE)
#   print(betabar)
#   cat(" ", fill = TRUE)
#   cat("A", fill = TRUE)
#   print(A)
#   cat(" ", fill = TRUE)
#   cat("dstarbar", fill = TRUE)
#   print(dstarbar)
#   cat(" ", fill = TRUE)
#   cat("Ad", fill = TRUE)
#   print(Ad)
#   cat(" ", fill = TRUE)
#   cat("MCMC parms: ", fill = TRUE)
#   cat("R= ", R, " keep= ", keep, "s= ", s, fill = TRUE)
#   cat(" ", fill = TRUE)
  betadraw = matrix(double(floor((R-burnin)/keep) * nvar), ncol = nvar)
  cutdraw = matrix(double(floor((R-burnin)/keep) * ncuts), ncol = ncuts)
  dstardraw = matrix(double(floor((R-burnin)/keep) * ndstar), ncol = ndstar)
  staydraw = array(0, dim = c((R-burnin)/keep))
  sigma = c(rep(1, nrow(X)))
  root = chol(chol2inv(chol((crossprod(X, X) + A))))
  Abetabar = crossprod(A, betabar)
  rootdi = chol(chol2inv(chol(Ad)))
  betahat = chol2inv(chol(crossprod(X, X))) %*% crossprod(X,y)
  dstarini = c(cumsum(c(rep(0.1, ndstar))))
  dstarout = optim(dstarini, lldstar, method = "BFGS", hessian = T, 
                   control = list(fnscale = -1, maxit = 500, reltol = 1e-06, 
                                  trace = 0), mu = X %*% betahat, y = y)
  inc.root = chol(chol2inv(chol((-dstarout$hessian + Ad))))
  olddstar = c(rep(0, ndstar))
  beta = betahat
  cutoffs = dstartoc(olddstar)
  oldll = lldstar(olddstar, y, mu = X %*% betahat)
#   itime = proc.time()[3]
#   cat("MCMC Iteration (est time to end - min) ", fill = TRUE)
#   fsh()
  Rnew = R-burnin
  withProgress(message = 'Making calculations', value = 0, {
  for (rep in 1:Rnew) { #ITER
    z = rtrun(X %*% beta, sigma = sigma, a = cutoffs[y], 
              b = cutoffs[y + 1])
    beta = breg1(root, X, z, Abetabar)
    metropout = dstarRwMetrop(y, X %*% beta, olddstar, s, 
                              inc.root, dstarbar, oldll, rootdi)
    olddstar = metropout$dstardraw
    oldll = metropout$oldll
    cutoffs = dstartoc(olddstar)
    stay = metropout$stay
#     if (rep%%100 == 0) {
#       ctime = proc.time()[3]
#       timetoend = ((ctime - itime)/rep) * (R - rep)
#       cat(" ", rep, " (", round(timetoend/60, 1), ")", 
#           fill = TRUE)
#       fsh()
#     }
    if (rep%%keep == 0) {
      mkeep = rep/keep
      cutdraw[mkeep, ] = cutoffs
      dstardraw[mkeep, ] = olddstar
      betadraw[mkeep, ] = beta
      staydraw[mkeep] = stay
    }
    incProgress(1/rep, detail = paste('Doing iteration', rep))
  }
  })
#  accept = 1 - sum(staydraw)/(Rnew/keep)
#   ctime = proc.time()[3]
#   cat("  Total Time Elapsed: ", round((ctime - itime)/60, 2), 
#       "\n")
  cutdraw = as.matrix(cutdraw[, 3:k])
  attributes(cutdraw)$class = "mcmc"
  attributes(betadraw)$class = "mcmc"
#  attributes(dstardraw)$class = "bayesm.mat"
  attributes(cutdraw)$mcpar = c(1, R, keep, burnin)
  attributes(betadraw)$mcpar = c(1, R, keep, burnin)
#  attributes(dstardraw)$mcpar = c(1, R, keep)
#   return(list(cutdraw = cutdraw, betadraw = betadraw, dstardraw = dstardraw, 
#               accept = accept))
return(list(betadraw = betadraw, cutdraw = cutdraw))
}


# setwd("C:/andres_ramirez/Andres_EAFIT/Bayesian Econometrics/UserInterface/BEsmartV21/Data")
# Data1<- read.table(file="SimOrderedmodel.csv",header=TRUE,sep=",")
# Data<-list(y=as.vector(Data1[,1]),X=as.matrix(Data1[,-1]))
# Bmean<-c(0,0,0,0)
# Bvar<- 0.01*diag(4)
# Prior<- list(Bmean,Bvar)
# MCMC<- list(R=1000,keep=20,burnin=1)
# 
# res<-Oprobit(Data,Prior,MCMC)
# res





