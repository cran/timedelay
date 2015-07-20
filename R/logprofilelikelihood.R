### log likelihood function of all the model parameters
logpostDelta <- function(delta, data, theta, c, log, unif) {

  time <- data[, 1]
  leng.time <- length(time)

  if (delta < unif[1] | delta > unif[2]) {

    -Inf

  } else {

    lcA <- data[, 2]
    se.lcA <- data[, 3]
    lcB <- data[, 4]
    se.lcB <- data[, 5]

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
    time.d <- time - delta
    time.temp <- c(time, time.d)
    ord <- order(time.temp)
    time.comb <- time.temp[ord]
    leng.time.comb <- length(time.comb)
  
    # indicator taking on 1 for X(t - delta) and 0 for X(t)
    ind <- c(rep(0, leng.time), rep(1, leng.time))[ord]

    lc.temp <- c(lcA, lcB - c)
    lc.comb <- lc.temp[ord]
    se.lc.temp <- c(se.lcA, se.lcB)
    se.lc.comb <- se.lc.temp[ord]
   
    if (all(se.lc.comb == 0)) {
      # Kelly et. al (2009)
      # x.star.i, i = 1, 2, ..., 2n
      x.star.i <- lc.comb - mu

      # omega.i, i = 1, 2, ..., 2n to be saved
      omega.i <- rep(NA, length(time.comb))

      # x.hat.i, i = 1, 2, ..., 2n to be saved
      x.hat.i <- rep(NA, length(time.comb))

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
      for (k in 2 : length(time.comb)) {
        x.hat.i[k] <- a.i[k - 1] * (x.hat.i[k - 1] +
                                    omega.i[k - 1] / (se.lc.comb[k - 1]^2 + omega.i[k - 1]) * 
                                    (x.star.i[k - 1] - x.hat.i[k - 1]))                  
      } 

      # log-likelihood
      sum(dnorm(x.star.i, mean = x.hat.i, sd = sqrt(omega.i + se.lc.comb^2), log = TRUE))

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



### log profile likelihood
logprofilelikelihood <- function(Delta, initial, data, data.flux, delta.uniform.range) {

  optim_delta <- function(th) {
    theta <- th[1 : 3]
    c <- th[4]
    logpostDelta(Delta, data, theta, c, log = data.flux, unif = delta.uniform.range) 
  }

  res.temp <- optim(initial, optim_delta, control = list(fnscale = -1), 
                    hessian = FALSE, method = "L-BFGS-B", lower = c(-30, 0, 0, -60), upper = c(30, Inf, Inf, 60))
  res.temp$value

}
