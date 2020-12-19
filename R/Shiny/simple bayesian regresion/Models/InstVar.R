InstVar<- function (Data, Prior, Mcmc) 
{
  pandterm = function(message) {
    stop(message, call. = FALSE)
  }
  if (missing(Data)) {
    pandterm("Requires Data argument -- list of z,w,x,y")
  }
  if (is.null(Data$z)) {
    pandterm("Requires Data element z")
  }
  z = Data$z
  if (is.null(Data$w)) {
    pandterm("Requires Data element w")
  }
  w = Data$w
  if (is.null(Data$x)) {
    pandterm("Requires Data element x")
  }
  x = Data$x
  if (is.null(Data$y)) {
    pandterm("Requires Data element y")
  }
  y = Data$y
  if (!is.vector(x)) {
    pandterm("x must be a vector")
  }
  if (!is.vector(y)) {
    pandterm("y must be a vector")
  }
  n = length(y)
  if (!is.matrix(w)) {
    pandterm("w is not a matrix")
  }
  if (!is.matrix(z)) {
    pandterm("z is not a matrix")
  }
  dimd = ncol(z)
  dimg = ncol(w)
  if (n != length(x)) {
    pandterm("length(y) ne length(x)")
  }
  if (n != nrow(w)) {
    pandterm("length(y) ne nrow(w)")
  }
  if (n != nrow(z)) {
    pandterm("length(y) ne nrow(z)")
  }
  if (missing(Prior)) {
    md = c(rep(0, dimd))
    Ad = 0.01 * diag(dimd)
    mbg = c(rep(0, (1 + dimg)))
    Abg = 0.01 * diag((1 + dimg))
    nu = 3
    V = diag(2)
  }
  else {
    if (is.null(Prior$md)) {
      md = c(rep(0, dimd))
    }
    else {
      md = Prior$md
    }
    if (is.null(Prior$Ad)) {
      Ad = 0.01 * diag(dimd)
    }
    else {
      Ad = Prior$Ad
    }
    if (is.null(Prior$mbg)) {
      mbg = c(rep(0, (1 + dimg)))
    }
    else {
      mbg = Prior$mbg
    }
    if (is.null(Prior$Abg)) {
      Abg = 0.01 * diag((1 + dimg))
    }
    else {
      Abg = Prior$Abg
    }
    if (is.null(Prior$nu)) {
      nu = 3
    }
    else {
      nu = Prior$nu
    }
    if (is.null(Prior$V)) {
      V = nu * diag(2)
    }
    else {
      V = Prior$V
    }
  }
  if (ncol(Ad) != nrow(Ad) || ncol(Ad) != dimd || nrow(Ad) != 
        dimd) {
    pandterm(paste("bad dimensions for Ad", dim(Ad)))
  }
  if (length(md) != dimd) {
    pandterm(paste("md wrong length, length= ", length(md)))
  }
  if (ncol(Abg) != nrow(Abg) || ncol(Abg) != (1 + dimg) || 
        nrow(Abg) != (1 + dimg)) {
    pandterm(paste("bad dimensions for Abg", dim(Abg)))
  }
  if (length(mbg) != (1 + dimg)) {
    pandterm(paste("mbg wrong length, length= ", length(mbg)))
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

  Rnew = R-burnin
  deltadraw = matrix(double(floor(Rnew/keep) * dimd), ncol = dimd)
  betadraw = rep(0, floor(Rnew/keep))
  gammadraw = matrix(double(floor(Rnew/keep) * dimg), ncol = dimg)
  Sigmadraw = matrix(double(floor(Rnew/keep) * 4), ncol = 4)
  Sigma = diag(2)
  delta = c(rep(0.1, dimd))
  xtd = matrix(nrow = 2 * n, ncol = dimd)
  ind = seq(1, (2 * n - 1), by = 2)
  zvec = as.vector(t(z))
  withProgress(message = 'Making calculations', value = 0, {
  for (rep in 1:Rnew) {
    e1 = as.vector(x - z %*% delta)
    ee2 = (Sigma[1, 2]/Sigma[1, 1]) * e1
    sig = sqrt(Sigma[2, 2] - (Sigma[1, 2]^2/Sigma[1, 1]))
    yt = (y - ee2)/sig
    xt = cbind(x, w)/sig
    bg = breg(yt, xt, mbg, Abg)
    beta = bg[1]
    gamma = bg[2:length(bg)]
    C = matrix(c(1, beta, 0, 1), nrow = 2)
    B = C %*% Sigma %*% t(C)
    L = t(chol(B))
    Li = backsolve(L, diag(2), upper.tri = FALSE)
    u = as.vector((y - w %*% gamma))
    yt = as.vector(Li %*% rbind(x, u))
    z2 = rbind(zvec, beta * zvec)
    z2 = Li %*% z2
    zt1 = z2[1, ]
    zt2 = z2[2, ]
    dim(zt1) = c(dimd, n)
    zt1 = t(zt1)
    dim(zt2) = c(dimd, n)
    zt2 = t(zt2)
    xtd[ind, ] = zt1
    xtd[-ind, ] = zt2
    delta = breg(yt, xtd, md, Ad)
    Res = cbind(x - z %*% delta, y - beta * x - w %*% gamma)
    S = crossprod(Res)
    Sigma = rwishart(nu + n, chol2inv(chol(V + S)))$IW

    if (rep%%keep == 0) {
      mkeep = rep/keep
      deltadraw[mkeep, ] = delta
      betadraw[mkeep] = beta
      gammadraw[mkeep, ] = gamma
      Sigmadraw[mkeep, ] = Sigma
    }
    incProgress(1/rep, detail = paste('Doing iteration', rep))
  }
  })

  attributes(deltadraw)$class = c("mcmc")
  attributes(deltadraw)$mcpar = c(1, R, keep, burnin)
  attributes(betadraw)$class = c("mcmc")
  attributes(betadraw)$mcpar = c(1, R, keep, burnin)
  attributes(gammadraw)$class = c("mcmc")
  attributes(gammadraw)$mcpar = c(1, R, keep, burnin)
  attributes(Sigmadraw)$class = c("mcmc")
  attributes(Sigmadraw)$mcpar = c(1, R, keep, burnin)
  return(list(deltadraw = deltadraw, betadraw = betadraw, gammadraw = gammadraw, 
              Sigmadraw = Sigmadraw))
}
