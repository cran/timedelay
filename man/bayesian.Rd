\name{bayesian}			
\alias{bayesian}
\title{Estimating the time delay via the Bayesian method}
\description{
  \code{bayesian} produces posterior samples of all the model parameters one of which is the time delay. This function  has options for three MCMC techniques, ancillarity-sufficiency interweaving strategy (ASIS), adaptive MCMC, and tempered transition, to improve the convergence rate of the MCMC and to handle the multimodality of the time delay.
}

\usage{bayesian(data, data.flux, theta.ini, 
   delta.ini, delta.uniform.range, delta.proposal.scale, 
   tau.proposal.scale, tau.prior.shape, tau.prior.scale,
   sigma.prior.shape, sigma.prior.scale, asis = TRUE,
   adaptive.freqeuncy, adaptive.delta =  TRUE, adaptive.delta.factor,
   adaptive.delta.acceptance.rate.upper.bound, 
   adaptive.delta.acceptance.rate.lower.bound,
   adaptive.tau = TRUE, adaptive.tau.factor,
   adaptive.tau.acceptance.rate.upper.bound,
   adaptive.tau.acceptance.rate.lower.bound,
   tempered.transition = FALSE, number.rungs, temperature.base,
   sample.size, warmingup.size)}

\arguments{
  \item{data}{
    A (\emph{n} by 1) matrix; the first column has n observation times, the second column has n flux (or magnitude) values of light A, the third column has n measurement errors of light A, the fourth column has n flux (or magnitude) values of light B, and the fifth column has n measurement errors of light B.
  }
  \item{data.flux}{
    "True" if data are recorded on flux scale or "FALSE" if data are on magnitude scale.  
  }
  \item{theta.ini}{
    Initial values of theta = (mu, sigma, tau) for MCMC.
  }
  \item{delta.ini}{
    Initial values of the time delay for MCMC.
  }
  \item{delta.uniform.range}{
    The range of the Uniform prior distribution for the time delay. The feasible entire support is c(min(simple[, 1]) - max(simple[, 1]), max(simple[, 1]) - min(simple[, 1])).
  }
  \item{delta.proposal.scale}{
    The proposal scale of the Metropolis step for the time delay.
  }
  \item{tau.proposal.scale}{
    The proposal scale of the Metropolis-Hastings step for tau.
  }
  \item{tau.prior.shape}{
    The shape parameter of the Inverse-Gamma hyper-prior distribution for tau. 
  }
  \item{tau.prior.scale}{
    The scale parameter of the Inverse-Gamma hyper-prior distribution for tau. 
  }
  \item{sigma.prior.shape}{
    The shape parameter of the Inverse-Gamma hyper-prior distribution for sigma^2.
  }
  \item{sigma.prior.scale}{
    The scale parameter of the Inverse-Gamma hyper-prior distribution for sigma^2. If no prior information is available, we recommend using 2 * 10^(-7).
  }
  \item{asis}{
    (Optional) "TRUE" if we use the ancillarity-sufficiency interweaving strategy (ASIS) for c (always recommended). Default is "TRUE".
  }
  \item{adaptive.freqeuncy}{
    (If "adaptive.delta = TRUE" or "adaptive.tau = TRUE") The adaptive MCMC is applied for every specified frequency. If it is specified as 500, the adaptive MCMC is applied to every 500th iterstion.
  }
  \item{adaptive.delta}{
    (Optional) "TRUE" if we use the adaptive MCMC for the time delay. Default is "TRUE".
  }
  \item{adaptive.delta.factor}{
    (If "adaptive.delta = TRUE") The factor, exp(adaptive.delta.factor) or exp(-adaptive.delta.factor), multiplied to the proposal scale of the time delay for adaptive MCMC.
  }
  \item{adaptive.delta.acceptance.rate.upper.bound}{
    (If "adaptive.delta = TRUE") We multiply exp(adaptive.delta.factor) to the proposal scale of the time delay if the acceptance rate of the time delay is greater than this number.
  }
  \item{adaptive.delta.acceptance.rate.lower.bound}{
    (If "adaptive.delta = TRUE") We multiply exp(-adaptive.delta.factor) to the proposal scale of the time delay if the acceptance rate of the time delay is smaller than this number.
  }
  \item{adaptive.tau}{
    (Optional) "TRUE" if we use the adaptive MCMC for tau. Default is "TRUE".
  }
  \item{adaptive.tau.factor}{
    (If "adaptive.tau = TRUE") The factor, exp(adaptive.tau.factor) or exp(-adaptive.tau.factor), multiplied to the proposal scale of tau for adaptive MCMC.
  }
  \item{adaptive.tau.acceptance.rate.upper.bound}{
    (If "adaptive.tau = TRUE") We multiply exp(adaptive.tau.factor) to the proposal scale of tau if the acceptance rate of tau is greater than this number.
  }
  \item{adaptive.tau.acceptance.rate.lower.bound}{
    (If "adaptive.tau = TRUE") We multiply exp(-adaptive.tau.factor) to the proposal scale of tau if the acceptance rate of tau is smaller than this number.
  }
  \item{tempered.transition}{
    (Optional) "TRUE" if we use the tempered transition. We recommend not using the adaptive MCMC for the time delay together with the tempered transition, i.e., "adaptive.delta = FALSE" if "tempered.transition = TRUE". Default is "FALSE".
  }
  \item{number.rungs}{
    (If "tempered.transition = TRUE") The number of rungs in the temperature ladder.
  }
  \item{temperature.base}{
    (If "tempered.transition = TRUE") The temerature of the i-th rung is (temperature.base)^i.
  }
  \item{sample.size}{
    The number of the posterior samples of each model parameter.
  }
  \item{warmingup.size}{
    The number of burn-ins for MCMC.
  }

}

\details{
  The function \code{bayesian} produces posterior samples of the model parameters one of which is the time delay.
}

\value{
  The outcome of \code{bayesian} comprises of:
  \item{delta}{Posterior samples of the time delay}
  \item{c}{Posterior samples of the magnitude offset, c}
  \item{mu}{Posterior samples of the mean parameter of the O-U process, mu}
  \item{sigma}{Posterior samples of the short term variability of the O-U process, sigma}
  \item{tau}{Posterior samples of the mean reversion time of the O-U process, tau}
  \item{tau.accept.rate}{The acceptance rate of the MCMC for tau}
  \item{delta.accept.rate}{The acceptance rate of the MCMC for the time delay}
}

\references{
    Hyungsuk Tak, Kaisey Mandel, David A. van Dyk, Vinay L. Kashyap, Xiao-Li Meng, and Aneta Siemiginowska (in progress). Bayesian and Profile Likelihood Approaches to Time Delay Estimation for Stochastic Time Series of Gravitationally Lensed Quasars
}

\author{
  Hyungsuk Tak
}
\examples{

  # Loading datasets
  data(simple)
  head(simple)

  ###############################################
  # Time delay estimation via Bayesian approach #
  ###############################################

  output = bayesian(dat = simple, data.flux = FALSE, theta.ini = c(0, 0.03, 100), 
                    delta.ini = 50, delta.uniform.range = c(0, 100), 
                    delta.proposal.scale = 0.3, 
                    tau.proposal.scale = 2, tau.prior.shape = 1, tau.prior.scale = 1,
                    sigma.prior.shape = 1, sigma.prior.scale = 2 * 10^(-7), asis = TRUE,
                    sample.size = 50, warmingup.size = 50)

  names(output)
  hist(output$delta, 50, xlab = expression(bold(Delta)), 
       main = expression(bold(paste("Histogram of posterior ", Delta))))
  plot(output$delta, ylab = expression(bold(Delta)), xlab = "Iteration",
       main = expression(bold(paste("Traceplot of posterior ", Delta))), type = "l")
  acf(output$delta, main = expression(bold(paste("ACF of posterior ", Delta))))


  ### tempered transition
  output = bayesian(dat = simple, data.flux = FALSE, theta.ini = c(0, 0.03, 100), 
                    delta.ini = 50, delta.uniform.range = c(0, 100), 
                    delta.proposal.scale = 5, 
                    tau.proposal.scale = 2, tau.prior.shape = 1, tau.prior.scale = 1,
                    sigma.prior.shape = 1, sigma.prior.scale = 2 * 10^(-7), asis = TRUE,
                    sample.size = 50, warmingup.size = 50, 
                    adaptive.delta = FALSE,
                    tempered.transition = TRUE, number.rungs = 2, temperature.base = 3)

  # We recommend using 10 rungs with the base of the temperature equal to 3.
  # We also recommend setting the "delta.proposal.scale" to 
  # the length of the feasible entire support / number of rungs

}


\keyword{methods}
