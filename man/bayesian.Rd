\name{bayesian}			
\alias{bayesian}
\title{Estimating the time delay via the Bayesian method}
\description{
  \code{bayesian} produces posterior samples of all the model parameters one of which is the time delay. This function has options for two MCMC techniques, ancillarity-sufficiency interweaving strategy (ASIS) and adaptive MCMC, to improve the convergence rate of the MCMC.
}

\usage{bayesian(data_lcA, data_lcB, data.flux, theta.ini, 
   delta.ini, delta.uniform.range, delta.proposal.scale, 
   tau.proposal.scale, tau.prior.shape, tau.prior.scale,
   sigma.prior.shape, sigma.prior.scale, asis = TRUE, micro,
   adaptive.freqeuncy, adaptive.delta =  TRUE, adaptive.delta.factor,
   adaptive.tau = TRUE, adaptive.tau.factor,
   sample.size, warmingup.size)}

\arguments{
  \item{data_lcA}{
    A (\emph{n} by 3) matrix; the first column has n observation times of light curve A, the second column has n flux (or magnitude) values of light curve A, the third column has n measurement errors of light curve A.
  }
  \item{data_lcB}{
    A (\emph{w} by 3) matrix; the first column has w observation times of light curve B, the second column has w flux (or magnitude) values of light curve B, the third column has w measurement errors of light curve B.
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
  \item{micro}{
    It determines the order of a polynomial regression model that accounts for the difference between microlensing trends. Default is 3. When zero is assigned, the Bayesian model fits a curve-shifted model.
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
  \item{adaptive.tau}{
    (Optional) "TRUE" if we use the adaptive MCMC for tau. Default is "TRUE".
  }
  \item{adaptive.tau.factor}{
    (If "adaptive.tau = TRUE") The factor, exp(adaptive.tau.factor) or exp(-adaptive.tau.factor), multiplied to the proposal scale of tau for adaptive MCMC.
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
  \item{beta}{Posterior samples of the polynomial regression coefficients, beta}
  \item{mu}{Posterior samples of the mean parameter of the O-U process, mu}
  \item{sigma}{Posterior samples of the short term variability of the O-U process, sigma}
  \item{tau}{Posterior samples of the mean reversion time of the O-U process, tau}
  \item{tau.accept.rate}{The acceptance rate of the MCMC for tau}
  \item{delta.accept.rate}{The acceptance rate of the MCMC for the time delay}
}

\references{
    Hyungsuk Tak, Kaisey Mandel, David A. van Dyk, Vinay L. Kashyap, Xiao-Li Meng, and Aneta Siemiginowska (2016+). "Bayesian Estimates of Astronomical Time Delays between Gravitationally Lensed Stochastic Light Curves," tentatively accepted in Annals of Applied Statistics (ArXiv 1602.01462).
}

\author{
  Hyungsuk Tak
}
\examples{

  # Loading datasets
  data(simple)

  # A typical quasar data set  
  head(simple)
  library(mnormt)

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
  output = bayesian(data_lcA = lcA, data_lcB = lcB, 
                    data.flux = FALSE, theta.ini = c(0, 0.01, 200), 
                    delta.ini = 50, delta.uniform.range = c(0, 100), 
                    delta.proposal.scale = 1, 
                    tau.proposal.scale = 3, 
                    tau.prior.shape = 1, tau.prior.scale = 1,
                    sigma.prior.shape = 1, sigma.prior.scale = 2 / 10^7, 
                    asis = TRUE, micro = 3,
                    sample.size = 100, warmingup.size = 50)

  names(output)

  # hist(output$delta, 20)
  # plot(output$delta, type = "l")
  # acf(output$delta)
  # output$delta.accept
  # output$tau.accept

}


\keyword{methods}
