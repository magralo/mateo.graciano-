XcreateMP<-function(nxs,nind,Data){
  pandterm = function(message) {
    stop(message, call. = FALSE)
  }
  if (missing(nxs)) 
    pandterm("requires number of regressors: include intercept if required")
  if (missing(nind)) 
    pandterm("requires number of units (individuals)")
  if (missing(Data)) 
    pandterm("requires dataset")
  if (nrow(Data)!=nind*2)
    pandterm("check dataset! number of units times number alternatives should be equal to dataset rows")
  p<-2
  XXDat<-array(0,c(p,1+nxs,nind))
  XX<-array(0,c(p,nxs*p,nind))
  YY<-array(0,c(p,1,nind))
  is<- seq(p,nind*p,p)
  cis<- seq(nxs,nxs*p+1,nxs)
  for(i in is){
    j<-which(i==is)
    XXDat[,,j]<-as.matrix(Data[c((i-(p-1)):i),-1])
    YY[,,j]<-XXDat[,1,j]
    for(l in 1:p){
      XX[l,((cis[l]-(nxs-1)):cis[l]),j]<-XXDat[l,-1,j]
    }
  }
  return(list(y=YY,X=XX))
}

BivProbit<-function(Data, Prior, Mcmc){
  pandterm = function(message) {
    print(Data)
    print(Prior)
    print(Mcmc)
    stop(message, call. = FALSE)
  }
  
  if (missing(Data)) {
    pandterm("Requires Data argument -- list of y, X")
  }
  if (is.null(Data$y)) {
    pandterm("Requires Data element y -- values of binary indicators")
  }
  y = Data$y
  if (is.null(Data$X)) {
    pandterm("Requires Data element X -- array of covariates")
  }
  X = Data$X
  p = dim(X)[1]
  k = dim(X)[2]
  N = dim(X)[3]
  ns = 1
  
  if (missing(Prior)) {
    betabar = rep(0, k)
    A = 100*diag(k)
    rhobar = 0
    v = 100
    tun = 0.25
  }
  else {
    if (is.null(Prior$betabar)) {
      betabar = rep(0, k)
    }
    else {
      betabar = Prior$betabar
    }
    if (is.null(Prior$A)) {
      A = 100*diag(k)
    }
    else {
      A = Prior$A
    }
    if (is.null(Prior$rhobar)) {
      rhobar = 0
    }
    else {
      rhobar = Prior$rhobar
    }
    if (is.null(Prior$v)) {
      v = 100
    }
    else {
      v = Prior$v
    }
    if (is.null(Prior$tun)) {
      tun = 0.25
    }
    else {
      tun = Prior$tun
    }
  }
  if (length(betabar) != k) 
    pandterm("length betabar ne k")
  if (sum(dim(A) == c(k, k)) != 2) 
    pandterm("A is of incorrect dimension")
  if (abs(rhobar)>1) 
    pandterm("invalid correlation coefficient. It should be between -1 and 1")
  if (tun<0 || tun>1) 
    pandterm("invalid tunning parameter in Random Walk Metropolis-Hastings algorithm. It should be between 0 and 1")
  if (missing(Mcmc)) 
    pandterm("Requires Mcmc argument -- at least R must be included")
  if (is.null(Mcmc$R)) {
    pandterm("Requires element R of Mcmc")
  }
  else {
    R = Mcmc$R
  }
  if (is.null(Mcmc$keep)) {
    pandterm("Requires element keep of Mcmc")
  }
  else {
    keep = Mcmc$keep
  }
  if (is.null(Mcmc$burnin)) {
    pandterm("Requires element burn-in of Mcmc")
  }
  else {
    burnin = Mcmc$burnin
  }
  
  betadraw = matrix(0,R,k)
  rhodraw = double(R)
  bet = c(rep(0, k))
  rho = 0
  Sig = matrix(c(1,rho,rho,1),2,2)
  mu<-matrix(0,p,N)
  yl<-matrix(0,p,N)
  Ainv<-solve(A)
  vsq<-v^0.5
  dfst<-3
  #Likelihood function
  log.L<-function(rhopt,yl,bet,X){
    N<-dim(X)[3]
    deter<- 1-rhopt^2
    Sig<- matrix(c(1,rhopt,rhopt,1),2,2)
    Sinv<-solve(Sig)
    logL<- NULL
    for(i in 1:N){
      logLi<-(-1/2)*log(deter)-(1/2)*t(yl[,i]-X[,,i]%*%bet)%*%Sinv%*%(yl[,i]-X[,,i]%*%bet)-(1/(2*v))*(rhopt-rhobar)^2
      logL<-c(logL,logLi)
    }
    logL<-sum(logL)
    return(logL)
  }
  withProgress(message = 'Making calculations', value = 0, {
  for (iter in 1:R) {
    for(h in 1:N){
      mu[,h]<-X[,,h]%*%bet
    }
    for (h in 1:N){
      if(y[1,1,h]==0){yl[1,h]<-rtruncnorm(1,-Inf,0,mu[1,h]+rho*(yl[2,h]-mu[2,h]),(1-rho^2)^0.5)} else{yl[1,h]<-rtruncnorm(1,0,Inf,mu[1,h]+rho*(yl[2,h]-mu[2,h]),(1-rho^2)^0.5)}
      if(y[2,1,h]==0){yl[2,h]<-rtruncnorm(1,-Inf,0,mu[2,h]+rho*(yl[1,h]-mu[1,h]),(1-rho^2)^0.5)} else{yl[2,h]<-rtruncnorm(1,0,Inf,mu[2,h]+rho*(yl[1,h]-mu[1,h]),(1-rho^2)^0.5)}
    }
    Sinv<-solve(Sig)
    AuxX<-array(0,c(k,k,N))
    Auxy<-array(0,c(k,1,N))
    for(h in 1:N){
      AuxX[,,h]<-t(X[,,h])%*%Sinv%*%X[,,h]
      Auxy[,,h]<-t(X[,,h])%*%Sinv%*%yl[,h]
    }
    AuxSumX<-apply(AuxX, MARGIN=c(1, 2), sum)
    AuxSumy<-apply(Auxy, MARGIN=c(1, 2), sum)
    B1<-solve(AuxSumX+Ainv)
    betmean<-B1%*%(AuxSumy+Ainv%*%betabar)
    bet<-mvrnorm(1,betmean,B1)
    repeat
    {
      rhoc<-rho+tun*rnorm(1)
      if (abs(rhoc)<1)
        break
    }
    alpha<-min(1,exp(log.L(rho=rhoc,yl=yl,bet=bet,X=X)-log.L(rho=rho,yl=yl,bet=bet,X=X)))
    u<-runif(1)
    if(u<=alpha){rho<-rhoc} else {rho<-rho}
    Sig<- matrix(c(1,rho,rho,1),2,2)
    betadraw[iter,]<-bet
    rhodraw[iter]<-rho
    incProgress(1/iter, detail = paste('Doing iteration', iter))
  }
  })
  Rnew = R-burnin
  betadrawMed<-betadraw[-c(1:burnin),]
  rhodrawMed<-rhodraw[-c(1:burnin)]
  betadrawNew=matrix(0,nrow=(floor(Rnew/keep)),ncol=k)
  rhodrawNew=double(floor(Rnew/keep))

  for (iter in 1:Rnew)
  {if (iter%%keep == 0) {
    mkeep = iter/keep
    betadrawNew[mkeep, ] = betadrawMed[iter,]
    rhodrawNew[mkeep] = rhodrawMed[iter]
  }}
  return(list(bet=mcmc(betadrawNew),rho=mcmc(rhodrawNew)))
}

