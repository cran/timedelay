
### log likelihood function of all the model parameters
logpostDelta <- function(delta, data_lcA, data_lcB, theta, c, log, unif, micro) {

  time1 <- data_lcA[, 1]
  time2 <- data_lcB[, 1]
  leng.time1 <- length(time1)
  leng.time2 <- length(time2)

  if (delta < unif[1] | delta > unif[2]) {

    -Inf

  } else if (theta[1] < -30 | theta[1] > 30) {

    -Inf

  } else {

    lcA <- data_lcA[, 2]
    se.lcA <- data_lcA[, 3]
    lcB <- data_lcB[, 2]
    se.lcB <- data_lcB[, 3]

    if (log == TRUE) {
      # transform into magnitude scale 
      se.lcA <- se.lcA * 2.5 / lcA / log(10)
      lcA <- -2.5 * log(lcA, base = 10)
      se.lcB <- se.lcB * 2.5 / lcB / log(10)
      lcB <- -2.5 * log(lcB, base = 10)
    }

    mu <- theta[1]
    sigma <- theta[2]
    tau <- theta[3]

    # sorting time given delta
    time.d <- time2 - delta
    time.temp <- c(time1, time.d)
    ord <- order(time.temp)
    time.comb <- time.temp[ord]
    leng.time.comb <- length(time.comb)

    # microlensing  
    if (micro == 0) {
      mat.temp <- matrix(c(rep(1, leng.time2)), ncol = 1)
    } else if (micro == 1) {
      mat.temp <- matrix(c(rep(1, leng.time2), time.d), ncol = 2)
    } else if (micro == 2) {
      mat.temp <- matrix(c(rep(1, leng.time2), time.d, time.d^2), ncol = 3)
    } else if (micro == 3) {
      mat.temp <- matrix(c(rep(1, leng.time2), time.d, time.d^2, time.d^3), 
                         ncol = 4)
    }

    c.pred <- mat.temp %*% c


    lc.temp <- c(lcA, lcB - c.pred)
    lc.comb <- lc.temp[ord]
    se.lc.temp <- c(se.lcA, se.lcB)
    se.lc.comb <- se.lc.temp[ord]
   
    if (all(se.lc.comb == 0)) {
      # Kelly et. al (2009)
      # x.star.i, i = 1, 2, ..., 2n
      x.star.i <- lc.comb - mu

      # omega.i, i = 1, 2, ..., 2n to be saved
      omega.i <- rep(NA, leng.time.comb)

      # x.hat.i, i = 1, 2, ..., 2n to be saved
      x.hat.i <- rep(NA, leng.time.comb)

      # a.i, i = 2, ..., 2n
      a.i <- exp( -diff(time.comb) / tau)

      # omega.i, i = 1, 2, ..., 2n
      omega.i[1] <- tau * sigma^2 / 2

      for (k in 2 : leng.time.comb) {
        omega.i[k] <- omega.i[1] * (1 - a.i[k - 1]^2) +
                      a.i[k - 1]^2 * omega.i[k - 1] * se.lc.comb[k - 1]^2 / 
                      (se.lc.comb[k - 1]^2 + omega.i[k - 1])
      }  

      # x.hat.i, i = 1, 2, ..., 2n
      x.hat.i[1] <- 0
      for (k in 2 : leng.time.comb) {
        x.hat.i[k] <- a.i[k - 1] * (x.hat.i[k - 1] +
                          omega.i[k - 1] / (se.lc.comb[k - 1]^2 + omega.i[k - 1]) * 
                          (x.star.i[k - 1] - x.hat.i[k - 1]))                  
      } 

      # log-likelihood
      sum(dnorm(x.star.i, mean = x.hat.i, sd = sqrt(omega.i + se.lc.comb^2), 
                log = TRUE))

    } else {

      # x.star.i, i = 1, 2, ..., 2n
      x <- lc.comb - mu

      # omega.i, i = 1, 2, ..., 2n to be saved
      B <- rep(NA, leng.time.comb)

      # x.hat.i, i = 1, 2, ..., 2n to be saved
      mu.i <- rep(NA, leng.time.comb)
      mu.star.i <- rep(NA, leng.time.comb)

      # a.i, i = 2, ..., 2n
      a.i <- exp( -diff(time.comb) / tau)

      # omega.i, i = 1, 2, ..., 2n
      var0 <- tau * sigma^2 / 2
      B[1] <- se.lc.comb[1]^2 / (se.lc.comb[1]^2 + var0)

      for (k in 2 : leng.time.comb) {
        B[k] <- se.lc.comb[k]^2 / ( se.lc.comb[k]^2 + 
                      a.i[k - 1]^2 * (1 - B[k - 1]) * se.lc.comb[k - 1]^2 +
                      var0 * (1 - a.i[k - 1]^2) ) 
      }  
  
      # x.hat.i, i = 1, 2, ..., 2n
      mu.i[1] <- (1 - B[1]) * x[1]
      for (k in 2 : leng.time.comb) {
        mu.i[k] <- (1 - B[k]) * x[k] + B[k] * a.i[k - 1] * mu.i[k - 1]
      }
  
      mu.star.i[1] <- 0
      mu.star.i[2 : leng.time.comb] <- a.i * mu.i[-leng.time.comb]

      var.star.i <- se.lc.comb^2 / B

      # log-likelihood
      sum(dnorm(x, mean = mu.star.i, sd = sqrt(var.star.i), log = TRUE))

    }  # end of if (all(se.lc.comb == 0)) 
  }  # end of if (delta < unif[1] | delta > unif[2])
}  # end of function



### entire log profile-likelihood curve
entirelogprofilelikelihood <- function(data_lcA, data_lcB, grid, 
                                       initial, data.flux, 
                                       delta.uniform.range, micro) {

  res.save <- rep(NA, length(grid))

  for (i in 1 : length(grid)) {

    delta.temp <- grid[i]

    optim_delta <- function(th) {
      mu <- th[1]
      sigma <- exp(th[2])
      tau <- exp(th[3])
      c <- th[4 : (micro + 4)]
      logpostDelta(delta = delta.temp, data_lcA = data_lcA, data_lcB = data_lcB,
                   theta = c(mu, sigma, tau), c = c, log = data.flux, 
                   unif = delta.uniform.range, micro = micro) 
    }

    res.temp <- optim(initial, optim_delta, control = list(fnscale = -1), 
                      hessian = FALSE, method = "BFGS")

    res.save[i] <- res.temp$value
    initial <- res.temp$par 

  }

  res.save
  
}
