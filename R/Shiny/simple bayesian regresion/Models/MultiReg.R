MultiReg<-function(Data, m, k, Prior, Mcmc){
  
  pandterm = function(message) {
    print(Data)
    print(m)
    print(k)
    print(Prior)
    print(Mcmc)
    stop(message, call. = FALSE)
  }
  
  if (missing(Data)) {
    pandterm("Requires Data argument -- list of Y and X")
  }
  if (missing(m)) {
    pandterm("Requires number of endogeneous varaibles")
  }
  if (missing(k)) {
    pandterm("Requires number of exogeneous varaibles")
  }
  
  if ((m+k) != ncol(Data)) {
    pandterm("Check data set! Dimensions do not agree")
  }
  
  Y = as.matrix(Data[,1:m])
  X = as.matrix(Data[,(m+1):(m+k)])
  
  if (missing(Prior)) {
    betabar=rep(0,k*m)
    Bbar=matrix(betabar,ncol=m)
    A=diag(rep(.01,k))
    nu=m-1
    V=nu*diag(m)
  }
  else {
    if (is.null(Prior$Bbar)) {
      Bbar = matrix(betabar,ncol=m)
    }
    else {
      Bbar = Prior$Bbar
    }
    if (is.null(Prior$A)) {
      A = diag(rep(.01,k))
    }
    else {
      A = Prior$A
    }
    if (is.null(Prior$nu)) {
      nu=m-1
    }
    else {
      nu = Prior$nu
    }
    if (is.null(Prior$V)) {
      V=nu*diag(m)
    }
    else {
      V = Prior$V
    }
  }
  
  if (missing(Mcmc)) {
    pandterm("Requires Mcmc argument")
  }
  else {
    if (is.null(Mcmc$R)) {
      pandterm("Requires Mcmc element R")
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
  
  
  betadraw=matrix(double(R*k*m),ncol=k*m)
  Sigmadraw=matrix(double(R*m*m),ncol=m*m)
  
  withProgress(message = 'Making calculations', value = 0, {
  for (rep in 1:R)
  {out=rmultireg(Y,X,Bbar,A,nu,V)
  betadraw[rep,]=out$B
  Sigmadraw[rep,]=out$Sigma
  incProgress(1/rep, detail = paste('Doing iteration', rep))
  }
  })
  
  Rnew = R-burnin
  betadrawMed<-betadraw[-c(1:burnin),]
  SigmadrawMed<-Sigmadraw[-c(1:burnin),]
  betadrawNew=matrix(0,nrow=(floor(Rnew/keep)),ncol=k*m)
  SigmadrawNew=matrix(0,nrow=(floor(Rnew/keep)),ncol=m*m)
  
  for (rep in 1:Rnew)
  {if (rep%%keep == 0) {
    mkeep = rep/keep
    betadrawNew[mkeep, ] = betadrawMed[rep,]
    SigmadrawNew[mkeep,] = SigmadrawMed[rep,]
  }}
  return(list(bet=mcmc(betadrawNew),Sig=mcmc(SigmadrawNew)))
}
