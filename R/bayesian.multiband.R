
Log.post.theta <- function(K, od, theta) {
  
  # theta: list, sigma, tau each is a K-vec, rho: K * K matrix
 
  sigma2 <- theta[[1]] # vector
  tau <- theta[[2]]  # actually is inverse tau!!
  rho <- theta[[3]] # K * K matrix: (1, rho, rho,1)
  
  
  if(K > 1){
    rho = matrix(rho, nrow = K, ncol =K)
    diag(rho) = 1 
  }else rho = matrix(rho)
  
  # sorting time given delta
  time.comb <- od$time.comb
  filter_ind <- od$filter_ind
  var.lc.comb  <- od$var
  lc.comb <- od$lc  # observations
  leng.time = length(time.comb)
  
  Sigma <- list() # T: K*K, sigma t|t-1
  Sigma.star <- list() # T: K*K, sigma t|t
  
  Sigma[[1]] = t(rho * sqrt(sigma2)) * sqrt(sigma2)  
  
  temp = matrix(0, K, K)
  for(i in 1:K)
  {
    for(j in 1:K){
      temp[i, j] = tau[i] + tau[j]
    }
  }
  Sigma[[1]]  =  Sigma[[1]] /temp # Q
  
  
  
  # x.hat.i, i = 1, 2, ..., 2n to be saved
  mu.i <- matrix(0, leng.time, K) # T * K, mu t|t-1
  mu.star.i <- matrix(0, leng.time, K) # mu t|t
  
  # a.i, i = 2, ..., 2n
  delta.t <- diff(time.comb)
  a.i <- exp( -delta.t %*% t(tau)) # (T-1) * K
  
  for (t in 1 : leng.time) {
    if(t>1) Sigma[[t]] <-  t(a.i[t - 1, ] *Sigma.star[[t-1]]) * a.i[t - 1, ] + Sigma[[1]] * (1 - exp(-temp*delta.t[t-1]))
    k =  ceiling(filter_ind[t]/2)
    Sigma.star[[t]] <- Sigma[[t]] - Sigma[[t]][, k] %*% t(Sigma[[t]][,k]) /(Sigma[[t]][k, k] + var.lc.comb[t])
    
    if(t>1) mu.i[t,] <- a.i[t - 1, ] * mu.star.i[t-1, ]
    mu.star.i[t,] <- mu.i[t,] + Sigma[[t]][, k] /(Sigma[[t]][k, k] + var.lc.comb[t]) * (lc.comb[t] - mu.i[t,k])
  }  
  
  # log-likelihood
  mus <- sapply(1:leng.time, function(t) mu.i[t, ceiling(filter_ind[t]/2)])
  vars <- sapply(1:leng.time, function(t) Sigma[[t]][ceiling(filter_ind[t]/2), ceiling(filter_ind[t]/2)] + var.lc.comb[t])
  if(any(vars<0))
  {
    print(vars[vars<0])
    print(time.comb[which(vars<0)])
  }
  return(sum(dnorm(lc.comb, mean = mus, sd = sqrt(vars), log = TRUE)))
}


postTheta.multi <- function (K, od, previous.theta, 
                       tau.thresh, sigma2.thresh, rho.thresh, 
                       tau.jump,sigma2.jump, rho.proposal.scale,
                       tau.prior.a, tau.prior.b, sigma.prior.a, 
                       sigma.prior.b, rho.prior.a, rho.prior.b) {
  
  # updating sigma
  sigma2 =  previous.theta[[1]]
  sigma2.p <- exp(log(sigma2) + sigma2.jump)
  for(k in 1:K)
  {
    theta = previous.theta
    theta[[1]][k] = sigma2.p[k]
    l.metrop <- Log.post.theta(K, od, theta) - Log.post.theta(K, od, previous.theta)
    l.metrop <- l.metrop - sigma.prior.b/sigma2.p[k] - (sigma.prior.a + 1) * log(sigma2.p[k]) 
                        + sigma.prior.b/sigma2[k] + (sigma.prior.a + 1) * log(sigma2[k])
    
    l.hastings <- log(sigma2.p[k]) - log(sigma2[k])
    
    # Accept-reject
    if (l.metrop + l.hastings > sigma2.thresh[k]) {
      previous.theta[[1]][k] <- sigma2.p[k]
    }
  }
  
  tau = previous.theta[[2]]
  # updating tau
  tau.p <- exp(log(tau) + tau.jump)
  for(k in 1:K)
  {
    theta = previous.theta
    theta[[2]][k] = tau.p[k]
    l.metrop <- Log.post.theta(K, od, theta) - Log.post.theta(K, od, previous.theta)
    l.metrop <- l.metrop - tau.p[k] * tau.prior.b + (tau.prior.a -1)*log(tau.p[k]) + tau[k] *tau.prior.b - (tau.prior.a -1)*log(tau[k])
    l.hastings <- log(tau.p[k]) - log(tau[k]) 
    
    # Accept-reject
    if (l.metrop + l.hastings > tau.thresh[k]) {
      previous.theta[[2]][k] <- tau.p[k]
    }
  }
  
  # updating rho
  if(K > 1)
  {
    rho =  previous.theta[[3]]
  inv.cdf <- runif(1, min = pnorm(-1, mean = rho, sd = rho.proposal.scale), 
                      max = pnorm(1, mean = rho, sd = rho.proposal.scale))
  rho.p <- qnorm(inv.cdf, mean = rho, sd = rho.proposal.scale)
    theta = previous.theta
    theta[[3]] = rho.p
    
    l.metrop <- Log.post.theta(K, od, theta) - Log.post.theta(K, od, previous.theta)
  l.hastings <- -log(pnorm((1 - rho.p) / rho.proposal.scale) - pnorm((-1 - rho.p) / rho.proposal.scale)) +
                 log(pnorm((1 - rho) / rho.proposal.scale) - pnorm((-1 - rho) / rho.proposal.scale))
    # Accept-reject
    if (l.metrop + l.hastings > rho.thresh) {
      previous.theta[[3]] <- rho.p
    }
  }
  
  return(previous.theta)  
  
}



#' posterior sampling of beta
postBeta <- function(od, X, K, theta, micro, hyperparam, timesca = 100, check = FALSE) # X: T *K
{
  time.d = od$time.comb
  filter_ind = od$filter_ind
  lc.comb = od$lc + od$c.pred  # original observations
  
  if(check) print(lc.comb[c(1,12,23,37,42)])
  
  c.pred = od$c.pred  # include mu
  
  rho =  theta[[3]] 
  sigma = sqrt(theta[[1]])
  tau = theta[[2]]
  
  leng.time = length(time.d)
  for(t in 1:leng.time)
    X[t, ceiling(filter_ind[t]/2)] = X[t, ceiling(filter_ind[t]/2)] + c.pred[t]
  
  delta.t <- diff(time.d)
  a.i <- exp( -delta.t %*% t(tau)) # (T-1) * K
  
  time.d.t = time.d/timesca
  if (micro == 0) {
    mat.temp <- matrix(c(rep(1, leng.time)), ncol = 1)
  } else if (micro == 1) {
    mat.temp <- matrix(c(rep(1, leng.time), time.d.t), ncol = 2)
  } else if (micro == 2) {
    mat.temp <- matrix(c(rep(1, leng.time), time.d.t, time.d.t^2), ncol = 3)
  } else if (micro == 3) {
    mat.temp <- matrix(c(rep(1, leng.time), time.d.t, time.d.t^2, time.d.t^3), ncol = 4)
  } else if (micro == 4) {
    mat.temp <- matrix(c(rep(1, leng.time), time.d.t, time.d.t^2, time.d.t^3, time.d.t^4), ncol = 5)
  } else if (micro == 5) {
    mat.temp <- matrix(c(rep(1, leng.time), time.d.t, time.d.t^2, time.d.t^3, time.d.t^4, time.d.t^5), ncol = 6)
  }
  
  Y = X[-1,] - a.i * X[-leng.time,] # T * K
  Y = rbind(X[1,] , Y)
  
  
  Q = t(rho * sigma) * sigma  
  
  tauM = matrix(0, K, K)
  for(i in 1:K)
  {
    for(j in 1:K){
      tauM[i, j] = tau[i] + tau[j]
    }
  }
  Q  =  Q /tauM   # Q
  
  micro = micro + 1
   
  Sigma_B = diag(rep(hyperparam[[2]], 2*K)) #original 1/hyperparam[[2]] , take the inverse after
  
  mB = matrix(rep(0, 2* micro * K),ncol=1)
  
  Sigma_B = cbind(Sigma_B, matrix(0, nrow(Sigma_B), K-1))
  Sigma_B = rbind(Sigma_B, matrix(0, K-1, ncol(Sigma_B)))
  
  Sigma_B[(2*micro*K + 1) : nrow(Sigma_B), (2*micro*K + 1) : nrow(Sigma_B)] = diag(rep(1, K-1), K-1,K-1)
  
  
  mB = rbind(mB, matrix(rep(0, K-1), K-1,1))
 
  for(t in 1:leng.time)
  {
    k = ceiling(filter_ind[t]/2)
    f = filter_ind[t]%%2
    W =  matrix(0, K, 2 * micro * K)
    I =  matrix(0, K, K - 1)
    
    if(f)
    {
      W1 = c(mat.temp[t, ] , rep(0,micro))
    }else{
      W1 = c(rep(0,micro), mat.temp[t, ] ) # 1 * 2m
    }
    
    
    if(t>1) 
    {
      if(delta.t[t-1] ==0) next
      
      k_1 = ceiling(filter_ind[t-1]/2)
      f_1 = filter_ind[t-1]%%2
      
      if(f_1)
      {
        W2 = c(-a.i[t-1, k_1] * mat.temp[t-1, ] , rep(0,micro))
      }else{
        W2 = c(rep(0,micro), -a.i[t-1, k_1] * mat.temp[t-1, ] ) # 1 * 2m
      }
      
      Qt = Q * (1 - exp(-tauM*delta.t[t-1]))
      a = a.i[t-1, k_1]
      

    }else{
      Qt = Q
      W2  = rep(0, 2 * micro)
      k_1 = k
      f_1 = f
      a = 0
    }
   
    inv_Qt = chol2inv(chol(Qt))
    
    range = (2 * (k-1) * micro +1) : (2 * k * micro)
    range1 = (2 * (k_1-1) * micro +1) : (2 * k_1 * micro)
    
    W[k, range] =  W[k, range] + W1
    W[k_1, range1] =  W[k_1, range1] + W2
    
    if(k > 1) I[k, k-1] = I[k, k-1] + 1
    if(k_1 > 1) I[k_1, k_1 -1] = I[k_1, k_1 -1] - a
    
    W = cbind(W, I)
    #print(min(eigen(Sigma_B, symmetric = T)$value))
    Sigma_B = Sigma_B + t(W)%*%inv_Qt%*%W
   
    mB= mB + t(Y[t, ] %*% inv_Qt %*% W)
    
  }
  
  #print(eigen(Sigma_B, symmetric = T)$value)
  Sigma_B =  chol2inv(chol(Sigma_B))
  mB = Sigma_B %*% mB
  
 
  B_comb = rmvnorm(1, mB, Sigma_B)
  
  if(any(is.na(B_comb)))
  {
    B_comb = t(mB)
  }
  
  
  beta = list()
  range=0
  for(k in 1:K)
  {
    beta[[k]] = matrix(B_comb[(range + 1) : (range + 2*micro)], ncol=2)
    range = range +  2*micro
  }
  
  
  c.pred <- sapply(1:leng.time,  function(i)  # T * 1 , will change if delta change!!!
  {
    k =  ceiling(filter_ind[i]/2)
    f = filter_ind[i] %% 2
    if(k>1) {
      return(mat.temp[i,] %*% beta[[k]][,2 - f] + B_comb[1, 2*micro*K +k -1])   # T * K = T * m , m * K
    }else{
      return(mat.temp[i,] %*% beta[[k]][,2 - f])
    }
   
  })
  
  
  lc.comb <- lc.comb - c.pred 
  
  mu = c(0, B_comb[1,-1:(-2*micro*K)])
  
  return(list("beta" = beta, "mu" = mu, "lc.comb" = lc.comb, "c.pred" = c.pred))
}

Log.post.delta <- function(delta, data_all, theta, c, c.mu, unif, micro) {

  # theta: list, sigma, tau each is a K-vec, rho: K * K matrix
  # c: list, beta, each is a m*2 matrix for m-degree poly and 2 curves
 
  K = length(data_all) # number of filters
  if (delta < unif[1] | delta > unif[2]) {
    -Inf
  } 
    sigma2 <- theta[[1]] # vector
    tau <- theta[[2]]  # actually is inverse tau!!
    rho <- theta[[3]] # K * K matrix: (1, rho, rho,1)
    
    
    od = orderData(data_all, delta, c, c.mu, micro = micro)

    time.comb = od$time.comb
    filter_ind = od$filter_ind
    var.lc.comb = od$var
    lc.comb = od$lc
    
    leng.time = length(time.comb)

    # omega.i, i = 1, 2, ..., 2n to be saved
    Sigma <- list() # T: K*K, sigma t|t-1
    Sigma.star <- list() # T: K*K, sigma t|t
    
    Sigma[[1]] = t(rho * sqrt(sigma2)) * sqrt(sigma2)  
    
    temp = matrix(0, K, K)
    for(i in 1:K)
    {
      for(j in 1:K){
        temp[i, j] = tau[i] + tau[j]
      }
    }
    Sigma[[1]]  =  Sigma[[1]] /temp # Q
    
   

    # x.hat.i, i = 1, 2, ..., 2n to be saved
    mu.i <- matrix(0, leng.time, K) # T * K, mu t|t-1
    mu.star.i <- matrix(0, leng.time, K) # mu t|t

    # a.i, i = 2, ..., 2n
    delta.t <- diff(time.comb)
    a.i <- exp( -delta.t %*% t(tau)) # (T-1) * K
    
    for (t in 1 : leng.time) {
      if(t>1) Sigma[[t]] <-  t(a.i[t - 1, ] *Sigma.star[[t-1]]) * a.i[t - 1, ] + Sigma[[1]] * (1 - exp(-temp*delta.t[t-1]))
      k =  ceiling(filter_ind[t]/2)
      Sigma.star[[t]] <- Sigma[[t]] - Sigma[[t]][, k] %*% t(Sigma[[t]][,k]) /(Sigma[[t]][k, k] + var.lc.comb[t])
      
      if(t>1) mu.i[t,] <- a.i[t - 1, ] * mu.star.i[t-1, ]
      mu.star.i[t,] <- mu.i[t,] + Sigma[[t]][, k] /(Sigma[[t]][k, k] + var.lc.comb[t]) * (lc.comb[t] - mu.i[t,k])
    }  

    # log-likelihood
    mus <- sapply(1:leng.time, function(t) mu.i[t, ceiling(filter_ind[t]/2)])
    vars <- sapply(1:leng.time, function(t) Sigma[[t]][ceiling(filter_ind[t]/2), ceiling(filter_ind[t]/2)] + var.lc.comb[t])
    if(any(vars<0))
    {
      print(vars[vars<0])
      print(time.comb[which(vars<0)])
    }
    return(sum(dnorm(lc.comb, mean = mus, sd = sqrt(vars), log = TRUE)))
}

postX.multi <- function(od, K, theta) { # input combined data remove polynomial
  
  # theta: list, sigma2, tau each is a K-vec, rho  K * K matrix
  time.comb = od$time.comb
  filter_ind = od$filter_ind
  var.lc.comb = od$var
  lc.comb = od$lc
  
  leng.time.comb = length(time.comb)
  

  sigma2 <- theta[[1]] # vector
  sigma = sqrt(sigma2)
  tau <- theta[[2]]  # actually is inverse tau!!
  rho <- theta[[3]] # K * K matrix: (1, rho, rho,1)
  

  Sigma <- list() # T: K*K, sigma t|t-1
  Sigma.star <- list() # T: K*K, sigma t|t
  
  Sigma[[1]] = t(rho * sigma) * sigma  
  
  temp = matrix(0, K, K)
  for(i in 1:K)
  {
    for(j in 1:K){
      temp[i, j] = tau[i] + tau[j]
    }
  }
  Sigma[[1]]  =  Sigma[[1]] /temp # Q
  
  
  
  # x.hat.i, i = 1, 2, ..., 2n to be saved
  mu.i <- matrix(0, leng.time.comb, K) # T * K, mu t|t-1
  mu.star.i <- matrix(0, leng.time.comb, K) # mu t|t
  
  # a.i, i = 2, ..., 2n
  delta.t <- diff(time.comb)
  a.i <- exp( -delta.t %*% t(tau)) # (T-1) * K
  
  for (t in 1 : leng.time.comb) {
    if(t>1) Sigma[[t]] <-  t(a.i[t - 1, ] *Sigma.star[[t-1]]) * a.i[t - 1, ] + Sigma[[1]] * (1 - exp(-temp*delta.t[t-1]))
    k =  ceiling(filter_ind[t]/2)
    Sigma.star[[t]] <- Sigma[[t]] - Sigma[[t]][, k] %*% t(Sigma[[t]][,k]) /(Sigma[[t]][k, k] + var.lc.comb[t])
    
    if(t>1) mu.i[t,] <- a.i[t - 1, ] * mu.star.i[t-1, ]
    mu.star.i[t,] <- mu.i[t,] + Sigma[[t]][, k] /(Sigma[[t]][k, k] + var.lc.comb[t]) * (lc.comb[t] - mu.i[t,k])
  }  
  
  
  # backward sampling
  X = matrix(0, leng.time.comb, K)
  mut = matrix(0, leng.time.comb, K)
  X[leng.time.comb,] <- rmvnorm(1, mu.star.i[leng.time.comb,], Sigma.star[[leng.time.comb]]) 
  mut[leng.time.comb,] <- mu.star.i[leng.time.comb,]
  for (t in  (leng.time.comb-1) :1) {
    mu_t = mu.star.i[t,] + Sigma.star[[t]] %*% ( chol2inv(chol(Sigma[[t+1]])) * a.i[t,])%*% (X[t+1, ] - a.i[t,]* mu.star.i[t,])
    sigma_t = Sigma.star[[t]] - Sigma.star[[t]] %*% (chol2inv(chol(Sigma[[t+1]]))* a.i[t,]) %*% (a.i[t,] * Sigma.star[[t]])
    
    mut[t, ] = mu.star.i[t,] + Sigma.star[[t]] %*% ( chol2inv(chol(Sigma[[t+1]])) * a.i[t,])%*% (mut[t+1, ] - a.i[t,]* mu.star.i[t,])
    if(delta.t[t] ==0) 
    {
      X[t, ] = X[t+1, ]
    }else{
      X[t, ]<- rmvnorm(1, mu_t, sigma_t)
    }
  }
  #X + c.pred 
  return(list('X' = X, 'mu' = mut))
}

#' @param c: microlensing polynomial coefficients
#' @param c.mu: mean diff of different bands
#' @return c.pred: polynomial part of the observations(including the means)
#' @return lc.comb: residual part
orderData <- function(data_all, delta, c, c.mu, timesc = 100, remove.poly = TRUE, micro = NULL)
{
  # sorting time given delta
  time.temp <- c()
  filter_ind <- c()
  var.lc.comb  <- c() 
  lc.comb <- c()  # observations
  
  k = 1
  for(data in data_all)
  {
    time <- data[, 1]
    time.d <- time - delta
    time.temp <- c(time.temp, time, time.d)
    filter_ind <-c(filter_ind, rep(k, length(time)), rep(k+1, length(time)))
    var.lc.comb  <- c(var.lc.comb , data[,3]^2, data[,5]^2)
    lc.comb <- c(lc.comb, data[,2], data[,4])
    k = k + 2
  }
  ord <- order(time.temp)
  time.comb <- time.temp[ord]
  filter_ind <- filter_ind[ord]
  var.lc.comb  <- var.lc.comb[ord]
  lc.comb <- lc.comb[ord]
  
  if(!remove.poly) return(list("time.comb" = time.comb, "filter_ind" = filter_ind, "var" = var.lc.comb, "lc" = lc.comb))
  
  leng.time <- length(time.comb)
  
  time.comb.t = time.comb/timesc
  if (micro == 0) {
    mat.temp <- matrix(c(rep(1, leng.time)), ncol = 1)
  } else if (micro == 1) {
    mat.temp <- matrix(c(rep(1, leng.time), time.comb.t), ncol = 2)
  } else if (micro == 2) {
    mat.temp <- matrix(c(rep(1, leng.time), time.comb.t, time.comb.t^2), ncol = 3)
  } else if (micro == 3) {
    mat.temp <- matrix(c(rep(1, leng.time), time.comb.t, time.comb.t^2, time.comb.t^3), ncol = 4)
  } 
  
  
  c.pred <- sapply(1:leng.time,  function(i)  # T * 1 , will change if delta change!!!
  {
    k =  ceiling(filter_ind[i]/2)
    f = filter_ind[i] %% 2
    return(mat.temp[i,] %*% c[[k]][,2 - f] + c.mu[k])   # T * K = T * m , m * K
    
  })
  
  
  lc.comb <- lc.comb - c.pred 
  
  return(list("time.comb" = time.comb, "filter_ind" = filter_ind, "var" = var.lc.comb, "lc" = lc.comb, "c.pred" = c.pred))
  
}







#' mcmc sampling
#' @param data_all: a list, each element is a data frame for one band, columns are time, value and std for different light curves (delayed). 
#' @param  theta.ini: initial values of the OU process parameters of multiple bands
#' @param delta.ini: initial value for time delay
###' @param kappa.ini: initial hyperparameter for the prior of microlensing polynomial coefficients
#' @param c.ini: initial coefficients of polynomials. A list, each element is a data frame for one band. Rows are coefficients of different degrees, columns are different light curves. Default: NULL
#' @param delta.uniform.range: prior of possible range of time delay
#' @param tau.*: paramters for prior, hyperprior and mcmc proposal distribution of 1/tau in OU process
#' @param sigma.*: paramters for prior, hyperprior and mcmc proposal distribution of sigma in OU process  
#' @param rho.*: paramters for prior and mcmc proposal distribution of rho in OU process  
###' @param kappa.*: paramters for prior and mcmc proposal distribution of kappa in the prior of microlensing polynomial coefficients
###' @param hyper.sigma2.*: hyperparameters for the microlensing polynomial coefficients
#' @param beta.prior.diag: prior for microlensing polynomial coefficients
#' @param micro: the degree of polynomial for microlensing
#' @param g: scale the time for fitting polynomial of microlensing
#' @param adaptive.*: parameters for adaptive MCMC
#' @param sample.size, warmingup,size: mcmc sample and burnin size


bayesian.multiband <- function(data.band1, data.band2, n.bands = 2, 
                               theta.ini = c(0.01, 0.01, 100, 100, 0.5), 
                               delta.ini, delta.uniform.range = c(-500, 500), delta.proposal.scale = 1, 
                               tau.proposal.scale = 1, tau.prior.shape = 1, tau.prior.scale = 1, 
                               sigma2.proposal.scale = 0.5, sigma2.prior.shape = 1, sigma2.prior.scale = 1e-7, 
                               rho.proposal.scale = 0.1,
                               beta.prior.diag = 10 * c(0.1, 0.01, 1e-3, 1e-5)^2,
                               micro = 3, timesc = 100, adaptive.frequency = 100,
                               adaptive.delta.factor = 0.1, adaptive.tau.factor = 0.1,
                               adaptive.sigma2.factor = 0.1, adaptive.rho.factor = 0.1,
                               sample.size, warmingup.size) {

  min.time <- min(data.band1[1, 1], data.band2[1, 1])
  data.band1[, 1] <- data.band1[, 1] - 3500
  data.band2[, 1] <- data.band2[, 1] - 3500
  K <- n.bands

  data_all <- vector(mode = "list", length = K)
  for (i in 1 : K) {
    eval(parse(text = paste0("data.temp <- data.band", i)))
    colnames(data.temp)[1] <- "V1"
    data_all[[i]] <- data.temp
  }

  beta.prior.diag <- beta.prior.diag[1 : (micro + 1)]

  total.sample.size <- sample.size + warmingup.size
  leng.time.comb = sum(sapply(data_all, nrow))
  
  sigma2.out <-  matrix(NA, total.sample.size, K)
  tau.out <- matrix(NA, total.sample.size, K)
  hyper.prior.out <- matrix(NA, total.sample.size, 4) # tau.shape, tau.scale, sigma2.shape, sigma2.scale
  delta.out <-  rep(NA, total.sample.size)
  rho.out <-  rep(NA, total.sample.size)
  c.out <- matrix(NA, nrow = total.sample.size, ncol = (micro + 1) * 2 * K)
  c.mu.out <- matrix(NA, nrow = total.sample.size, ncol = K)
  
  X.out <- matrix(NA, nrow = total.sample.size, ncol = leng.time.comb * K*2)  
  loglik.out <- rep(NA, total.sample.size)
  
  tau.accept <- matrix(0, total.sample.size, K)
  sigma2.accept <- matrix(0, total.sample.size, K)
  delta.accept <- rep(0, total.sample.size)
  rho.accept <- rep(0, total.sample.size)

  tau.jumps <- tau.proposal.scale * matrix(rnorm(total.sample.size * K), nrow = K, ncol =total.sample.size)
  tau.thresh <- matrix(-rexp(total.sample.size * K), ncol = K, nrow =total.sample.size)
  tau.proposal.scale.adapt <- rep(1,K)
  
  sigma2.jumps <- sigma2.proposal.scale * matrix(rnorm(total.sample.size * K), nrow = K, ncol =total.sample.size)
  sigma2.thresh <- matrix(-rexp(total.sample.size * K), ncol = K, nrow =total.sample.size) 
  sigma2.proposal.scale.adapt <- rep(1,K)

  delta.jumps <- delta.proposal.scale * rnorm(total.sample.size)
  delta.thresh <- -rexp(total.sample.size)
  delta.proposal.scale.adapt <- 1
  
  rho.thresh <- -rexp(total.sample.size)
  rho.proposal.scale.adapt <- rho.proposal.scale
  
  delta.t <- delta.ini  # delta ini
  sigma2.t <- theta.ini[1 : 2]^2  #sigma ini
  tau.t <- 1 / theta.ini[3 : 4]  #tau ini
  rho <- theta.ini[5]
  
  if(K > 1){
    rho.t = matrix(rho, nrow = K, ncol =K) 
    diag(rho.t) = 1 
  }else rho.t = matrix(rho)
  
  theta.prior.t = list(tau.prior.shape, tau.prior.scale, sigma2.prior.shape, sigma2.prior.scale)
  hyperparam.t = list(0, beta.prior.diag, rep(0, micro+1))
  
  c.t = list()
  c.mu.t = rep(0, K)
  i =1
  mu0 = mean(c(data_all[[1]][,2], data_all[[1]][,4]))  # mean of different light curves for the first band 
  for(dat in data_all)
  {
    c.mu.t[i] = mean(c(dat[,2], dat[,4])) - mu0 # different in mean relatively to the first band 
    ti1 <- dat[, 1]/timesc
    ti2 <- ti1^2
    ti3 <- ti1^3

    if (micro == 0) {
      ss1 <- lm(dat[, 2] - c.mu.t[i]  ~ 1) 
    } else if (micro == 1) {
      ss1 <- lm(dat[, 2] - c.mu.t[i]  ~ ti1) 
    } else if (micro == 2) {
      ss1 <- lm(dat[, 2] - c.mu.t[i]  ~ ti1 + ti2) 
    } else if (micro == 3) {
      ss1 <- lm(dat[, 2] - c.mu.t[i]  ~ ti1 + ti2 + ti3) 
    }
      
    ti1 <- (dat[, 1] - delta.ini)/timesc
    ti2 <- ti1^2
    ti3 <- ti1^3
    if (micro == 0) {
      ss2 <- lm(dat[, 4] - c.mu.t[i] ~ 1) 
    } else if (micro == 1) {
      ss2 <- lm(dat[, 4] - c.mu.t[i] ~ ti1) 
    } else if (micro == 2) {
      ss2 <- lm(dat[, 4] - c.mu.t[i] ~ ti1 + ti2) 
    } else if (micro == 3) {
      ss2 <- lm(dat[, 4] - c.mu.t[i] ~ ti1 + ti2 + ti3) 
    }

    if(micro > 3){
      c.t[[i]] <- cbind(c(ss1$coefficients, rep(0,micro - 2)), c(ss2$coefficients, rep(0,micro - 2)))
    }else{
      c.t[[i]] <- cbind(ss1$coefficients, ss2$coefficients)
    }
    i = i +1
  }
  
  print(paste("Starting time:", Sys.time()))
  od <- orderData(data_all, delta.t, c.t, c.mu.t, remove.poly = T, micro = micro) # remove the means of different bands/curves
  
  for (i in 1 : total.sample.size) {
  
    reslist <- postX.multi(od, K, theta = list(sigma2.t, tau.t, rho.t))
    X.t = reslist$X
    mu.t = reslist$mu
    
    tmp_list = postBeta(od, X.t, K, theta = list(sigma2.t, tau.t, rho.t), micro, hyperparam.t, timesca =timesc, check = FALSE) # X: T *K
      
    c.t = tmp_list$beta
    od$lc.comb = tmp_list$lc.comb
    c.mu.t = tmp_list$mu
    
    # theta update
    tau.jump.adapt <- tau.proposal.scale.adapt * tau.jumps[,i]
    sigma2.jump.adapt <- sigma2.proposal.scale.adapt * sigma2.jumps[,i]
    
    
    if(K > 1) rho = rho.t[1,2]
    theta.update <- postTheta.multi(K, od, previous.theta = list(sigma2.t, tau.t, rho), 
                              tau.jump = tau.jump.adapt, tau.thresh = tau.thresh[i, ],
                              tau.prior.a = tau.prior.shape, tau.prior.b = tau.prior.scale, 
                              sigma2.jump = sigma2.jump.adapt, sigma2.thresh = sigma2.thresh[i, ],
                              sigma.prior.a =  sigma2.prior.shape, sigma.prior.b = sigma2.prior.scale,
                              rho.proposal.scale = rho.proposal.scale.adapt, rho.thresh = rho.thresh[i])
    
     if (K > 1){
       if( theta.update[[3]] != rho.t[1,2]) {
      rho.accept[i] <- 1
      rho.t = matrix(theta.update[[3]], nrow = K, ncol =K)
      diag(rho.t) = 1
      
       }
     }
    
    tau.accept[i, theta.update[[2]] != tau.t] <- 1
    sigma2.accept[i, theta.update[[1]] != sigma2.t] <- 1

    # delta update
    delta.p <- delta.t + delta.proposal.scale.adapt * delta.jumps[i]
    prev_lik <- Log.post.delta(delta.t, data_all, list(sigma2.t, tau.t, rho.t), c.t, c.mu.t, 
                               unif = delta.uniform.range, micro)
    cur_lik <- Log.post.delta(delta.p, data_all, list(sigma2.t, tau.t, rho.t), c.t, c.mu.t, 
                              unif = delta.uniform.range, micro)
    l.metrop <- cur_lik - prev_lik
    
    if (l.metrop > delta.thresh[i]) { 
      delta.t <- delta.p 
      delta.accept[i] <- 1
      od <- orderData(data_all, delta.t, c.t, c.mu.t, remove.poly = TRUE, micro = micro)
      
      loglik.out[i] <- cur_lik
    }else loglik.out[i] <- prev_lik
    
   
    sigma2.t <- sigma2.out[i,] <- theta.update[[1]]
    tau.t <- tau.out[i,] <- theta.update[[2]]
    rho.out[i] <- theta.update[[3]]
    
    X.out[i,] <- mu.t
    c.out[i, ]<- unlist(c.t)
    c.mu.out[i, ] <- c.mu.t
    delta.out[i] <- delta.t
        
    if (i %% adaptive.frequency == 0) {
      if(mean(delta.accept[i - (adaptive.frequency - 1) : i]) > 0.35) {
        scale.adj <- exp(min(adaptive.delta.factor, 1 / sqrt(i / adaptive.frequency)))
      } else if (mean(delta.accept[i - (adaptive.frequency - 1) : i]) < 0.35) {
        scale.adj <- exp(-min(adaptive.delta.factor, 1 / sqrt(i / adaptive.frequency)))
      } else {
        scale.adj <- 1
      }
      delta.proposal.scale.adapt <- delta.proposal.scale.adapt * scale.adj
        
    }
    
    if (i %% adaptive.frequency == 0) {
      if(mean(rho.accept[i - (adaptive.frequency - 1) : i]) > 0.35) {
        scale.adj <- exp(min(adaptive.delta.factor, 1 / sqrt(i / adaptive.frequency)))
      } else if (mean(rho.accept[i - (adaptive.frequency - 1) : i]) < 0.35) {
        scale.adj <- exp(-min(adaptive.delta.factor, 1 / sqrt(i / adaptive.frequency)))
      } else {
        scale.adj <- 1
      }
      rho.proposal.scale.adapt <- rho.proposal.scale.adapt * scale.adj
        
    }
    
    if (i %% adaptive.frequency == 0) {
      scale.adj <- exp(min(adaptive.tau.factor, 1 / sqrt(i / adaptive.frequency)))
      
      ind = which(colMeans(tau.accept[(i - (adaptive.frequency - 1)) : i, , drop=FALSE]) > 0.35)
      tau.proposal.scale.adapt[ind] <- tau.proposal.scale.adapt[ind] * scale.adj
        
      ind = which(colMeans(tau.accept[(i - (adaptive.frequency - 1)) : i, , drop=FALSE]) < 0.35)
      tau.proposal.scale.adapt[ind] <- tau.proposal.scale.adapt[ind] / scale.adj
        
    }
    
    if (i %% adaptive.frequency == 0) {
      scale.adj <- exp(min(adaptive.sigma2.factor, 1 / sqrt(i / adaptive.frequency)))
      
      ind = which(colMeans(sigma2.accept[(i - (adaptive.frequency - 1)) : i, , drop=FALSE]) > 0.35)
      sigma2.proposal.scale.adapt[ind] <- sigma2.proposal.scale.adapt[ind] * scale.adj
        
      ind = which(colMeans(sigma2.accept[(i - (adaptive.frequency - 1)) : i, , drop=FALSE]) < 0.35)
      sigma2.proposal.scale.adapt[ind] <- sigma2.proposal.scale.adapt[ind] / scale.adj
        
    }
    
  }
  
  print(paste("Ending time:", Sys.time()))
  
  out <- list(delta = delta.out[-c(1 : warmingup.size)],
              beta = c.out[-c(1 : warmingup.size), ],
              rho = rho.out[-c(1 : warmingup.size)], 
              sigma = sqrt(sigma2.out[-c(1 : warmingup.size),]), 
              tau = 1 / tau.out[-c(1 : warmingup.size),],
              tau.accept.rate = colMeans(tau.accept),
              sigma.accept.rate = colMeans(sigma2.accept),
              delta.accept.rate = mean(delta.accept),
              rho.accept.rate = mean(rho.accept))
  
  out
    

}
