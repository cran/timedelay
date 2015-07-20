\name{logprofilelikelihood}			
\alias{logprofilelikelihood}
\title{Calculating the profilel likelihood function}
\description{
  \code{logprofilelikelihood} calculates the value of the profile likelihood function given a specific value of the time delay. 
}

\usage{logprofilelikelihood(Delta, initial, data, data.flux, delta.uniform.range)}

\arguments{
  \item{Delta}{
    A value of the time delay at which the profile likelihood value is calculated.
  }
  \item{initial}{
    The initial values of the other model parameters (mu, sigma, tau, c)
  }
  \item{data}{
    The data set   
  }
  \item{data.flux}{
    "True" if data are recorded on flux scale or "FALSE" if data are on magnitude scale.
  }
  \item{delta.uniform.range}{
    The range of the Uniform prior distribution for the time delay. The feasible entire support is c(min(simple[, 1]) - max(simple[, 1]), max(simple[, 1]) - min(simple[, 1])).
  }

}

\details{
  The function \code{logprofilelikelihood} can be used to find the maximum profile likelihood estimate and the inverse of the negative Hessian at the mode. See the example below.
}

\value{
  The outcome of \code{logprofilelikelihood} is the value of the log profile likelihood function at the given time delay.
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

  ################################################
  # Time delay estimation via profile likelihood #
  ################################################

  ###### Individual profile likelihood value
  theta.ini = c(0, 0.03, 100)
  c.ini <- mean(simple[, 4]) - mean(simple[, 2])
  delta.uniform.range <- c(0, 100)

  ###  Running the following codes takes more time than CRAN policy
  ###  Please type the following two lines without "#"  to run the function
  #  logprofilelikelihood(Delta = 40, initial = c(theta.ini, c.ini), data = simple, 
  #                       data.flux = FALSE, delta.uniform.range = delta.uniform.range)

  ###### Maximum profile likelihood estimate and 
  ###### the inverse of the negative Hessian value at the mode
  
  # prof.lik <- function(delta) {
  #   logprofilelikelihood(delta, initial = c(theta.ini, c.ini), data = simple, 
  #                        data.flux = FALSE, delta.uniform.range = delta.uniform.range)
  # }
  # outcome <- optim(50, prof.lik, control = list(fnscale = -1), method = "BFGS", hessian = TRUE)

  # maximum profile likelihood estimate:
  # outcome$par
  # inverse of the negative Hessian value at the mode:
  # 1 / sqrt(-outcome$hessian)


}


\keyword{methods}
