\name{bayesian}			
\alias{bayesian}
\title{Estimating the time delay via the Bayesian method}
\description{
  \code{bayesian} produces posterior samples of all the model parameters including the time delay. This function has options for three Markov chain Monte Carlo (MCMC) techniques; (i) ancillarity-sufficiency interweaving strategy (ASIS) and (ii) adaptive MCMC to improve the convergence rate, and (iii) repelling-attracting Metropolis for multimodality. 
}

\usage{bayesian(data.lcA, data.lcB, data.flux, theta.ini, 
   delta.ini, delta.uniform.range, delta.proposal.scale, 
   tau.proposal.scale, tau.prior.shape = 1, tau.prior.scale = 1,
   sigma.prior.shape = 1, sigma.prior.scale = 1e-7, 
   asis = TRUE, micro, multimodality = FALSE, 
   adaptive.frequency, adaptive.delta =  TRUE, adaptive.delta.factor,
   adaptive.tau = TRUE, adaptive.tau.factor,
   sample.size, warmingup.size)}

\arguments{
  \item{data.lcA}{
    A (\eqn{n_1} by 3) matrix for image A (light curve A); the first column has \eqn{n_1} observation times, the second column contains \eqn{n_1} flux (or magnitude) values, the third column includes \eqn{n_1} measurement errors.
  }
  \item{data.lcB}{
    A (\eqn{n_2} by 3) matrix for image B (light curve B); the first column has \eqn{n_2} observation times, the second column contains \eqn{n_2} flux (or magnitude) values, the third column includes \eqn{n_2} measurement errors.
  }
  \item{data.flux}{
    "True" if data are recorded on flux scale or "FALSE" if data are on magnitude scale. Simply speaking, if your observed time series can take on negative values, then data.flux = FALSE.
  }
  \item{theta.ini}{
    Initial values of the OU parameters, (\eqn{\mu}, \eqn{\sigma}, \eqn{\tau}).
  }
  \item{delta.ini}{
    Initial value of the time delay, \eqn{\Delta}.
  }
  \item{delta.uniform.range}{
    The range of the Uniform prior distribution for the time delay.
  }
  \item{delta.proposal.scale}{
    The proposal scale of the Metropolis step for the time delay.
  }
  \item{tau.proposal.scale}{
    The proposal scale of the Metropolis-Hastings step for \eqn{\tau}.
  }
  \item{tau.prior.shape}{
    The shape parameter of the Inverse-Gamma hyper-prior distribution for \eqn{\tau}. Default is 1.
  }
  \item{tau.prior.scale}{
    The scale parameter of the Inverse-Gamma hyper-prior distribution for \eqn{\tau}. Default is 1.
  }
  \item{sigma.prior.shape}{
    The shape parameter of the Inverse-Gamma hyper-prior distribution for \eqn{\sigma^2}.  Default is 1.
  }
  \item{sigma.prior.scale}{
    The scale parameter of the Inverse-Gamma hyper-prior distribution for \eqn{\sigma^2}. Default is 1e-7.
  }
  \item{asis}{
    (Optional) "TRUE" if we use the ancillarity-sufficiency interweaving strategy (ASIS) for c (always recommended). Default is "TRUE".
  }
  \item{micro}{
    It determines the order of a polynomial regression model that accounts for the difference between microlensing trends. Default is 3. When zero is assigned, the Bayesian model fits a curve-shifted model.
  }
  \item{multimodality}{
    (Optional) "TRUE" if we use a repelling-attracting Metropolis algorithm for sampling Delta. When it is "TRUE" it is recommended that adaptive.delta = FALSE so that the adaptive MCMC is not used in the presence of multimodality. Default is "FALSE". 
  }
  \item{adaptive.frequency}{
    (If "adaptive.delta = TRUE" or "adaptive.tau = TRUE") The adaptive MCMC is applied for every specified frequency. If it is specified as 500, the adaptive MCMC is applied to every 500th iterstion.
  }
  \item{adaptive.delta}{
    (Optional) "TRUE" if we use the adaptive MCMC for the time delay. Default is "TRUE".
  }
  \item{adaptive.delta.factor}{
    (If "adaptive.delta = TRUE") The factor, exp(\eqn{\pm}adaptive.delta.factor), multiplied by the proposal scale of the time delay for adaptive MCMC.
  }
  \item{adaptive.tau}{
    (Optional) "TRUE" if we use the adaptive MCMC for tau. Default is "TRUE".
  }
  \item{adaptive.tau.factor}{
    (If "adaptive.tau = TRUE") The factor, exp(\eqn{\pm}adaptive.tau.factor), multiplied by the proposal scale of tau for adaptive MCMC.
  }
  \item{sample.size}{
    The number of the posterior samples of each model parameter.
  }
  \item{warmingup.size}{
    The number of burn-in posterior samples.
  }

}

\details{
  The function \code{bayesian} produces posterior samples of the model parameters one of which is the time delay. Please note that when astronomical time series data are loaded on R by \code{read.table}, \code{read.csv}, etc., some decimal places of the the observation times are automatically rounded because R's default is to load seven effective digits. For example, R will load the observation time 51075.412789 as 51075.41. This default will produce many ties in observation times even though there is actually no tie in observation times. To prevent this, please type "\code{options(digits = 11)}" before loading the data if the observation times are in seven effective digits.
}

\value{
  The outcome of \code{bayesian} comprises of:
  \item{delta}{Posterior samples of the time delay}
  \item{X}{Posterior samples of the latent magnitudes}
  \item{beta}{Posterior samples of the polynomial regression coefficients, beta}
  \item{mu}{Posterior samples of the mean parameter of the O-U process, mu}
  \item{sigma}{Posterior samples of the short term variability of the O-U process, sigma}
  \item{tau}{Posterior samples of the mean reversion time of the O-U process, tau}
  \item{tau.accept.rate}{The acceptance rate of the MCMC for tau}
  \item{delta.accept.rate}{The acceptance rate of the MCMC for the time delay}
}

\references{
    Hyungsuk Tak, Kaisey Mandel, David A. van Dyk, Vinay L. Kashyap, Xiao-Li Meng, and Aneta Siemiginowska (2017). "Bayesian Estimates of Astronomical Time Delays between Gravitationally Lensed Stochastic Light Curves," The Annals of Applied Statistics, 11 (3), 1309-1348.
    Hyungsuk Tak, Xiao-Li Meng, and David A. van Dyk (2018), "A Repelling-Attracting Metropolis Algorithm for Multimodality", Journal of Computational and Graphical Statistics, 27 (3), 479-490.
    Zhirui Hu and Hyungsuk Tak (2020+), "Modeling Stochastic Variability in  Multi-Band Time Series Data," arXiv:2005.08049.
}

\author{
  Hyungsuk Tak
}
\examples{

  # Loading datasets
  data(simple)

  # A typical quasar data set  
  head(simple)

  # Subset (data for image A) of the typical quasar data set
  lcA <- simple[, 1 : 3]

  # Another subset (data for image B) of the typical quasar data set
  # The observation times for image B are not necessarily the same as those for image A
  lcB <- simple[, c(1, 4, 5)]

  # The two subsets do not need to have the same number of observations
  # For example, here we add one more observation time for image B
  lcB <- rbind(lcB, c(290, 1.86, 0.006))

  dim(lcA)
  dim(lcB)

  ###############################################
  # Time delay estimation via Bayesian approach #
  ###############################################

  # Cubic microlensing model (m = 3)
  \donttest{
  output = bayesian(data.lcA = lcA, data.lcB = lcB, 
                    data.flux = FALSE, theta.ini = c(0, 0.01, 200), 
                    delta.ini = 50, delta.uniform.range = c(0, 100), 
                    delta.proposal.scale = 1, 
                    tau.proposal.scale = 3, 
                    tau.prior.shape = 1, tau.prior.scale = 1,
                    sigma.prior.shape = 1, sigma.prior.scale = 2 / 10^7, 
                    asis = TRUE, micro = 3,
                    sample.size = 10, warmingup.size = 10)
  names(output)
  }
  # hist(output$delta, 20)
  # plot(output$delta, type = "l")
  # acf(output$delta)
  # output$delta.accept
  # output$tau.accept

  # If multimodality exists, then please add two arguments:
  # multimodality = TRUE, adaptive.delta = FALSE
  # This is to prevent the Markov chain from adapting to a local mode.
  # If we know the distance between the most distant modes,
  # try the smallest value of delta.proposal.scale that still makes 
  # the Markov chain jump between those modes frequently.
  # For example, when the distance is d, then try 
  # delta.proposal.scale = d * 0.7, decreasing or increasing 0.7.

  # Graphical model checking 
  # beta.hat <- colMeans(output$beta)
  # delta.hat <- mean(output$delta)
  # time.x <- lcB[, 1] - delta.hat
  # time.covariate <- cbind(rep(1, length(time.x)), time.x, time.x^2, time.x^3)
  ##### This time.covariate is when micro = 3, a third-order polynomial regression.
  ##### With micro = 0, "time.covariate <- rep(1, length(time.x))" is used.
  # predicted <- time.covariate %*% beta.hat
  # plot(lcA[, 1], lcA[, 2], col = 2)     
  ##### Adjust the range of the x-axis by the argument of plot.
  ##### For example "xlim = c(-1, 2)" 
  # points(lcB[, 1] - delta.hat, lcB[, 2] - predicted, col = 4, pch = 1)
  # for (i in 1 : length(output$delta)) {
  #   temp.time <- c(lcA[, 1], lcB[, 1] - output$delta[i])
  #   points(sort(temp.time), output$X[i, ], 
  #          col = "gray", cex = 0.5, pch = 1, lwd = 0.01)
  # }
  # points(lcA[, 1], lcA[, 2], pch = 0, col = 2, lwd = 1)
  # points(lcB[, 1] - delta.hat, lcB[, 2] - predicted, 
  #        col = 4, pch = 1, lwd = 1)


}


\keyword{methods}
