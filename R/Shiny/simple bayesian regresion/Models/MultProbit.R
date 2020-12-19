MultinomialProbit<- function (Data, Prior, Mcmc) 
{
  pandterm = function(message) {
    stop(message, call. = FALSE)
  }
  if (missing(Data)) {
    pandterm("Requires Data argument -- list of p, y, X")
  }
  if (is.null(Data$p)) {
    pandterm("Requires Data element p -- number of alternatives")
  }
  p = Data$p
  if (is.null(Data$y)) {
    pandterm("Requires Data element y -- number of alternatives")
  }
  y = Data$y
  if (is.null(Data$X)) {
    pandterm("Requires Data element X -- matrix of covariates")
  }
  X = Data$X
  levely = as.numeric(levels(as.factor(y)))
  if (length(levely) != p) {
    pandterm(paste("y takes on ", length(levely), " values -- must be ", 
                   p))
  }
  bady = FALSE
  for (i in 1:p) {
    if (levely[i] != i) 
      bady = TRUE
  }
#   cat("Table of y values", fill = TRUE)
#   print(table(y))
  if (bady) {
    pandterm("Invalid y")
  }
  n = length(y)
  k = ncol(X)
  pm1 = p - 1
  if (nrow(X)/n != pm1) {
    pandterm(paste("X has ", nrow(X), " rows; must be = (p-1)n"))
  }
  if (missing(Prior)) {
    betabar = rep(0, k)
    A = 0.01 * diag(k)
    nu = pm1 + 3
    V = nu * diag(pm1)
  }
  else {
    if (is.null(Prior$betabar)) {
      betabar = rep(0, k)
    }
    else {
      betabar = Prior$betabar
    }
    if (is.null(Prior$A)) {
      A = 0.01 * diag(k)
    }
    else {
      A = Prior$A
    }
    if (is.null(Prior$nu)) {
      nu = pm1 + 3
    }
    else {
      nu = Prior$nu
    }
    if (is.null(Prior$V)) {
      V = nu * diag(pm1)
    }
    else {
      V = Prior$V
    }
  }
  if (length(betabar) != k) 
    pandterm("length betabar ne k")
  if (sum(dim(A) == c(k, k)) != 2) 
    pandterm("A is of incorrect dimension")
  if (nu < 1) 
    pandterm("invalid nu value")
  if (sum(dim(V) == c(pm1, pm1)) != 2) 
    pandterm(paste("V is of incorrect dimension",pm1))
  if (missing(Mcmc)) 
    pandterm("Requires Mcmc argument -- at least R must be included")
  if (is.null(Mcmc$R)) {
    pandterm("Requires element R of Mcmc")
  }
  else {
    R = Mcmc$R
  }
  if (is.null(Mcmc$beta0)) {
    beta0 = rep(0, k)
  }
  else {
    beta0 = Mcmc$beta0
  }
  if (is.null(Mcmc$sigma0)) {
    sigma0 = diag(pm1)
  }
  else {
    sigma0 = Mcmc$sigma0
  }
  if (length(beta0) != k) 
    pandterm("beta0 is not of length k")
  if (sum(dim(sigma0) == c(pm1, pm1)) != 2) 
    pandterm("sigma0 is of incorrect dimension")
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
#   cat(" ", fill = TRUE)
#   cat("Starting Gibbs Sampler for MNP", fill = TRUE)
#   cat("  ", n, " obs; ", p, " choice alternatives; ", k, " indep vars (including intercepts)", 
#       fill = TRUE)
#   cat("  ", R, " reps; keeping every ", keep, "th draw", fill = TRUE)
#   cat(" ", fill = TRUE)
#   cat("Table of y values", fill = TRUE)
#   print(table(y))
#   cat("Prior Parms:", fill = TRUE)
#   cat("betabar", fill = TRUE)
#   print(betabar)
#   cat("A", fill = TRUE)
#   print(A)
#   cat("nu", fill = TRUE)
#   print(nu)
#   cat("V", fill = TRUE)
#   print(V)
#   cat(" ", fill = TRUE)
#   cat("MCMC Parms:", fill = TRUE)
#   cat("R= ", R, fill = TRUE)
#   cat("initial beta= ", beta0, fill = TRUE)
#   cat("initial sigma= ", sigma0, fill = TRUE)
#   cat(" ", fill = TRUE)
  Rnew = R-burnin
  sigmadraw = matrix(double(floor(Rnew/keep) * pm1 * pm1), ncol = pm1 * 
                       pm1)
  betadraw = matrix(double(floor(Rnew/keep) * k), ncol = k)
  wnew = double(nrow(X))
  betanew = double(k)
  wold = c(rep(0, nrow(X)))
  betaold = beta0
  C = chol(solve(sigma0))
  Rcpp::sourceCpp("Models/draww.cpp")
  draww = function(w, X, y, beta, sigmai) {
    Xbeta = as.vector(X %*% beta)
    drawwc(w, Xbeta, y, sigmai)
  }
  drawwc = function(w, mu, y, sigi) {
    .C("draww", w = as.double(w), as.double(mu), as.double(sigi), 
       as.integer(length(y)), as.integer(ncol(sigi)), as.integer(y))$w
  }
  
 
  
#   itime = proc.time()[3]
#   cat("MCMC Iteration (est time to end - min) ", fill = TRUE)
  ## ITER
  withProgress(message = 'Making calculations', value = 0, {
  for (rep in 1:Rnew) {
    sigmai = crossprod(C)
    wnew = draww(wold, X, y, betaold, sigmai)
    zmat = matrix(cbind(wnew, X), nrow = pm1)
    zmat = C %*% zmat
    zmat = matrix(zmat, nrow = nrow(X))
    betanew = breg(zmat[, 1], zmat[, 2:(k + 1)], betabar, 
                   A)
    epsilon = matrix((wnew - X %*% betanew), nrow = pm1)
    S = crossprod(t(epsilon))
    W = rwishart(nu + n, chol2inv(chol(V + S)))
    C = W$C
#     if (rep%%100 == 0) {
#       ctime = proc.time()[3]
#       timetoend = ((ctime - itime)/rep) * (R + 1 - rep)
#       cat(" ", rep, " (", round(timetoend/60, 1), ")", 
#           fill = TRUE)
#       fsh()
#     }
    if (rep%%keep == 0) {
      mkeep = rep/keep
      betadraw[mkeep, ] = betanew
      sigmadraw[mkeep, ] = as.vector(W$IW)
    }
    wold = wnew
    betaold = betanew
    incProgress(1/rep, detail = paste('Doing iteration', rep))
  }
  })
#   ctime = proc.time()[3]
#   cat("  Total Time Elapsed: ", round((ctime - itime)/60, 2), 
#       "\n")
  attributes(betadraw)$class = c("mcmc")
  attributes(betadraw)$mcpar = c(1, R, keep, burnin)
  attributes(sigmadraw)$class = c("mcmc")
  attributes(sigmadraw)$mcpar = c(1, R, keep, burnin)
#  list(betadraw = betadraw/sigmadraw[,1]^0.5, sigmadraw = sigmadraw/sigmadraw[,1])
  list(betadraw = betadraw, sigmadraw = sigmadraw)
}



