\name{bayesian.multiband}			
\alias{bayesian.multiband}
\title{Estimating the time delay between doubly-lensed multi-band light curves in a Bayesian way}
\description{
\code{bayesian.multiband} produces posterior samples of all the model parameters including the time delay between doubly-lensed multi-band light curves.
}

\usage{bayesian.multiband(data.band1, data.band2, n.bands = 2,
  theta.ini = c(0.01, 0.01, 100, 100, 0.5),
  delta.ini, delta.uniform.range = c(-500, 500), delta.proposal.scale = 1,
  tau.proposal.scale = 1, tau.prior.shape = 1, tau.prior.scale = 1,
  sigma2.proposal.scale = 0.5, sigma2.prior.shape = 1, sigma2.prior.scale = 1e-7,
  rho.proposal.scale = 0.1, beta.prior.diag = 10 * c(0.1, 0.01, 1e-3, 1e-5)^2,
  micro = 3, timesc = 100, adaptive.frequency = 100,
  adaptive.delta.factor = 0.1, adaptive.tau.factor = 0.1,
  adaptive.sigma2.factor = 0.1, adaptive.rho.factor = 0.1,
  sample.size, warmingup.size)}

\arguments{
  \item{data.band1}{
    A (\eqn{n_1} by 5) matrix composed of the data for doubly-lensed images oberved in band 1; the first column contains \eqn{n_1} observation times, the second column contains \eqn{n_1} magnitudes for image A, the third column has \eqn{n_1} measurement error standard deviations for image A, the fourth column include \eqn{n_1} magnitudes for image B, and the fifth column contains \eqn{n_1} measurement error standard deviations for image B.
  }
  \item{data.band2}{
    A (\eqn{n_2} by 5) matrix composed of the data for doubly-lensed images oberved in band 2; the first column contains \eqn{n_2} observation times, the second column contains \eqn{n_2} magnitudes for image A, the third column has \eqn{n_2} measurement error standard deviations for image A, the fourth column include \eqn{n_2} magnitudes for image B, and the fifth column contains \eqn{n_2} measurement error standard deviations for image B.
  }
  \item{n.bands}{
    The number of bands used to obtain the data. For now, only two bands are allowed. The authors plan to modify the code to allow more bands in the near future.
  }
  \item{theta.ini}{
    A vector for initial values of the OU processes for two bands, i.e., (\eqn{\sigma_1}, \eqn{\sigma_2}, \eqn{\tau_1}, \eqn{\tau_2}, \eqn{\rho}), where \eqn{\sigma_1} is the short-term variability (standard deviation) for band 1, \eqn{\sigma_2} is the short-term variability (standard deviation) for band 2, \eqn{\tau_1} is the timescale for band 1, and \eqn{\tau_2} is the timescale for band 2, and \eqn{\rho} is the cross-correlation parameter between two bands. Defaults are (\eqn{0.01}, \eqn{0.01}, 100, 100, 0.5).
  }
  \item{delta.ini}{
    Initial values of the time delay.
  }
  \item{delta.uniform.range}{
    The range of the Uniform prior distribution for the time delay. Default range is set to \eqn{(-500, 500)}.
  }
  \item{delta.proposal.scale}{
    The proposal scale of the Metropolis step for the time delay. Default is 1.
  }
  \item{tau.proposal.scale}{
    The proposal scale of the Metropolis-Hastings step for \eqn{\tau_j} (\eqn{j=1, 2}) Default is 1.
  }
  \item{tau.prior.shape}{
    The shape parameter of the Inverse-Gamma hyper-prior distribution for \eqn{\tau_j}. Default is 1.
  }
  \item{tau.prior.scale}{
    The scale parameter of the Inverse-Gamma hyper-prior distribution for \eqn{\tau_j}. Default is 1.
  }
  \item{sigma2.proposal.scale}{
    The proposal scale of the Metropolis-Hastings step for \eqn{\sigma^2_j} (\eqn{j=1, 2}). Default is 0.5.
  }
  \item{sigma2.prior.shape}{
    The shape parameter of the Inverse-Gamma hyper-prior distribution for \eqn{\sigma^2_j}. Default is 1.
  }
  \item{sigma2.prior.scale}{
    The scale parameter of the Inverse-Gamma hyper-prior distribution for \eqn{\sigma^2_j}. If no prior information is available, we recommend using 1e-7 (default).
  }
  \item{rho.proposal.scale}{
    The proposal scale of the Metropolis-Hastings step for \eqn{\rho}. Default is 0.1.
  }
  \item{beta.prior.diag}{
    The diagonal elements of the covariance matrix in the multivariate Gaussian prior for \eqn{\beta} (polynomial regression coefficients for microlensing adjustment). If such information is not available, these are set to \code{10 * c(0.1)^2} if \code{micro = 0}, \code{10 * c(0.1, 0.01)^2} if \code{micro = 1}, \code{10 * c(0.1, 0.01, 1e-3)^2} if \code{micro = 2}, and \code{10 * c(0.1, 0.01, 1e-3, 1e-5)^2} if \code{micro = 3}.
  }
  \item{micro}{
    A non-negative integer less than or equal to 3. It determines the order of a polynomial regression model that accounts for the long-term trend of microlensing effect. Default is 3.
  }
  \item{timesc}{
    It scales the observation time for fitting polynomial of microlensing, i.e., time / timesc, so that the coefficients are not too small. Default is 100.
  }
  \item{adaptive.frequency}{
    The adaptive MCMC is applied for every specified frequency. If it is specified as 500, the adaptive MCMC is applied to every 500th iterstion. Default is 100.
  }
  \item{adaptive.delta.factor}{
    The factor, exp(\eqn{\pm}adaptive.delta.factor), multiplied by the proposal scale of the time delay for adaptive MCMC. Default is 0.1.
  }
  \item{adaptive.tau.factor}{
    The factor, exp(\eqn{\pm}adaptive.tau.factor), multiplied by the proposal scale of \eqn{\tau_j} for adaptive MCMC. Default is 0.1.
  }
  \item{adaptive.sigma2.factor}{
    The factor, exp(\eqn{\pm}adaptive.tau.factor), multiplied by the proposal scale of \eqn{\sigma^2_j} for adaptive MCMC. Default is 0.1.
  }
  \item{adaptive.rho.factor}{
    The factor, exp(\eqn{\pm}adaptive.tau.factor), multiplied by the proposal scale of \eqn{\rho} for adaptive MCMC. Default is 0.1.
  }
  \item{sample.size}{
    The number of the posterior samples of each model parameter.
  }
  \item{warmingup.size}{
    The number of burn-in samples for MCMC.
  }

}

\details{
The function \code{bayesian.multiband} produces posterior samples of the model parameters, where the time delay (delta) is of primary interest. For now, this function only supports doubly-lensed data observed in two bands. The authors plan to generalize this code to account for more than two bands and more than two lens.

Please note that when astronomical time series data are loaded on R by \code{read.table}, \code{read.csv}, etc., some decimal places of the the observation times are automatically rounded because R's default is to load seven effective digits. For example, R will load the observation time 51075.412789 as 51075.41. This default will produce many ties in observation times even though there is actually no tie in observation times. To prevent this, please type "\code{options(digits = 11)}" before loading the data if the observation times are in seven effective digits.
}

\value{
  The outcomes of \code{bayesian.multiband} are composed of:
  \item{delta}{A vector for \eqn{m} posterior samples of the time delay.}
  \item{beta}{An \eqn{m} by \eqn{k + 1} matrix containing posterior samples of the polynomial regression coefficients, where \eqn{m} is the size of the posterior sample, and \eqn{k} is the polynomial order for microlensing.}
  \item{rho}{A vector for \eqn{m} posterior samples of the cross-correlation parameter.}
  \item{sigma}{An \eqn{m} by \eqn{2} matrix containing posterior samples of the short-term variability (standard deviation) of the O-U process. The first column is composed ot the \eqn{m} posterior samples of \eqn{\sigma_1} in band 1, and the second column contains the \eqn{m} posterior samples of \eqn{\sigma_2} in band 2.}
  \item{tau}{An \eqn{m} by \eqn{2} matrix containing posterior samples of the timescale of the O-U process. The first column is composed ot the \eqn{m} posterior samples of \eqn{\tau_1} in band 1, and the second column contains the \eqn{m} posterior samples of \eqn{\tau_2} in band 2.}
  \item{tau.accept.rate}{The acceptance rate of the MCMC for \eqn{\tau_1} and \eqn{\tau_2}.}
  \item{sigma.accept.rate}{The acceptance rate of the MCMC for \eqn{\sigma_1} and \eqn{\sigma_2}.}
  \item{delta.accept.rate}{The acceptance rate of the MCMC for the time delay}
  \item{rho.accept.rate}{The acceptance rate of the MCMC for \eqn{\rho}.}
}

\references{
    Hyungsuk Tak, Kaisey Mandel, David A. van Dyk, Vinay L. Kashyap, Xiao-Li Meng, and Aneta Siemiginowska (2017). "Bayesian Estimates of Astronomical Time Delays between Gravitationally Lensed Stochastic Light Curves," The Annals of Applied Statistics, 11 (3), 1309-1348.
    Hyungsuk Tak, Xiao-Li Meng, and David A. van Dyk (2018), "A Repelling-Attracting Metropolis Algorithm for Multimodality", Journal of Computational and Graphical Statistics, 27 (3), 479-490.
}

\author{
  Zhirui Hu and Hyungsuk Tak
}
\examples{

  # Loading datasets 
  data(simple.band1)
  data(simple.band2)

  # Doubly-lensed quasar data set observed in two bands
  # Each data set contains doubly-lensed light curves observed in one band. 
  head(simple.band1)
  head(simple.band2)

  # The length of each data set (i.e., number of observation times)
  # do not need to be the same.

  dim(simple.band1)
  dim(simple.band2)

  ###############################################
  # Time delay estimation via Bayesian approach #
  ###############################################

  # Cubic microlensing model (m = 3)
  \donttest{
  output <- bayesian.multiband(data.band1 = simple.band1, 
              data.band2 = simple.band2, n.bands = 2, 
              theta.ini = c(0.01, 0.01, 100, 100, 0.5),
              delta.ini = 100, delta.uniform.range = c(-500, 500), 
              tau.proposal.scale = 1, tau.prior.shape = 1, tau.prior.scale = 1, 
              sigma2.proposal.scale = 0.5, sigma2.prior.shape = 1, sigma2.prior.scale = 1e-7, 
              rho.proposal.scale = 0.1, beta.prior.diag = 10 * c(0.1, 0.01, 1e-3, 1e-5)^2, 
              micro = 3, timesc = 100, adaptive.frequency = 100,
              adaptive.delta.factor = 0.1, adaptive.tau.factor = 0.1,
              adaptive.sigma2.factor = 0.1, adaptive.rho.factor = 0.1,
              sample.size = 100, warmingup.size = 100)
  names(output)
  }

  # hist(output$delta, 20)
  # plot(output$delta, type = "l")
  # acf(output$delta)
  # output$delta.accept.rate

}

\keyword{methods}
