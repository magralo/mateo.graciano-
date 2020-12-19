SUR<-function (Data, Prior, Mcmc) 
{
  pandterm = function(message) {
    print(Data)
    print(Prior)
    print(Mcmc)
    stop(message, call. = FALSE)
  }
  if (missing(Data)) {
    pandterm("Requires Data argument -- list of regdata")
  }
  if (is.null(Data$regdata)) {
    pandterm("Requires Data element regdata")
  }
  regdata = Data$regdata
  nreg = length(regdata)
  nobs = length(regdata[[1]]$y)
  nvar = 0
  indreg = double(nreg + 1)
  y = NULL
  for (reg in 1:nreg) {
    if (length(regdata[[reg]]$y) != nobs || nrow(regdata[[reg]]$X) != 
        nobs) {
      pandterm(paste("incorrect dimensions for regression", 
                     reg))
    }
    else {
      indreg[reg] = nvar + 1
      nvar = nvar + ncol(regdata[[reg]]$X)
      y = c(y, regdata[[reg]]$y)
    }
  }
  indreg[nreg + 1] = nvar + 1
  if (missing(Prior)) {
    betabar = c(rep(0, nvar))
    A = 0.01 * diag(nvar)
    nu = nreg-1
    V = nu * diag(nreg)
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
    if (is.null(Prior$nu)) {
      nu = nreg-1
    }
    else {
      nu = Prior$nu
    }
    if (is.null(Prior$V)) {
      V = nu * diag(nreg)
    }
    else {
      V = Prior$V
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
  }
  cat(" ", fill = TRUE)
  cat("Starting Gibbs Sampler for SUR Regression Model", fill = TRUE)
  cat("  with ", nreg, " regressions", fill = TRUE)
  cat("  and  ", nobs, " observations for each regression", 
      fill = TRUE)
  
  nk = ncol(regdata[[1]]$X)
  Xstar = regdata[[1]]$X
  Y = regdata[[1]]$y
  
  for (i in 2:nreg) {
    nk[i] = ncol(regdata[[i]]$X)
    Xstar = bdiag(Xstar, regdata[[i]]$X)
    Y = c(Y, regdata[[i]]$y)
  }
  k = sum(nk)
  tXstar<- t(Xstar)

  
  E = matrix(double(nobs * nreg), ncol = nreg)
  for (reg in 1:nreg) {
    E[, reg] = lm(y ~ . - 1, data = data.frame(y = regdata[[reg]]$y, 
                                               regdata[[reg]]$X))$residuals
  }
  Sigma = (crossprod(E) + diag(0.01, nreg))/nobs
  Sigmainv = chol2inv(chol(Sigma))
  Abetabar = A %*% betabar
  Vinv = solve(V)
  vNew<-nu+nobs
  
  betadraw=matrix(0,nrow=R,ncol=k)
  Sigmadraw=matrix(0,nrow=R,ncol=nreg*nreg)

  withProgress(message = 'Making calculations', value = 0, {
    for (rep in 1:R){
      B = solve(tXstar%*%kronecker(diag(nobs),Sigmainv)%*%Xstar + A)
      btilde = B%*%(tXstar%*%kronecker(diag(nobs),Sigmainv)%*%Y + Abetabar)
      beta = mvrnorm(n = 1, btilde, B)
      betadraw[rep,]=matrix(beta,1,k)
      res<-matrix(Y-Xstar%*%beta,nobs,nreg)
      ScaleIw<-solve(Vinv+t(res)%*%res)
      Sigmainv<-matrix(rWishart(1, vNew, ScaleIw),nreg,nreg)
      Sigmadraw[rep,]=c(solve(Sigmainv))
    incProgress(1/rep, detail = paste('Doing iteration', rep))
      }
  })
  
  Rnew = R-burnin
  betadrawMed<-betadraw[-c(1:burnin),]
  SigmadrawMed<-Sigmadraw[-c(1:burnin),]
  betadrawNew=matrix(0,nrow=(floor(Rnew/keep)),ncol=k)
  SigmadrawNew=matrix(0,nrow=(floor(Rnew/keep)),ncol=nreg*nreg)
  

  for (rep in 1:Rnew)
  {if (rep%%keep == 0) {
    mkeep = rep/keep
    betadrawNew[mkeep, ] = betadrawMed[rep,]
    SigmadrawNew[mkeep,] = SigmadrawMed[rep,]
  }}
  return(list(bet=mcmc(betadrawNew),Sig=mcmc(SigmadrawNew)))
}

  