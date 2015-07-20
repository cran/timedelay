\name{entirelogprofilelikelihood}			
\alias{entirelogprofilelikelihood}
\title{Calculating the entire profilel likelihood curve over the given grid values of the time delay}
\description{
  \code{entirelogprofilelikelihood} calculates the entire profilel likelihood curve over the given grid values of the time delay. 
}

\usage{entirelogprofilelikelihood(data, grid, initial, data.flux, delta.uniform.range)}

\arguments{
  \item{data}{
    The data set   
  }
  \item{grid}{
    A vector containing values of the time delay on which the profile likelihood values are calculated. We recommend using the grid interval equal to 0.1.
  }
  \item{initial}{
    The initial values of the other model parameters (mu, sigma, tau, c)
  }
  \item{data.flux}{
    "True" if data are recorded on flux scale or "FALSE" if data are on magnitude scale.
  }
  \item{delta.uniform.range}{
    The range of the Uniform prior distribution for the time delay. The feasible entire support is c(min(simple[, 1]) - max(simple[, 1]), max(simple[, 1]) - min(simple[, 1])).
  }

}

\details{
  The function \code{entirelogprofilelikelihood} is used to obtain the entire profile likelihood curve over the given grid values of the time delay. 
}

\value{
  The outcome of \code{entirelogprofilelikelihood} is the values of the log profile likelihood function over the given grid values of the time delay.
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


  ###### The entire profile likelihood values on the grid of values of the time delay.

  theta.ini = c(0, 0.03, 100)
  c.ini <- mean(simple[, 4]) - mean(simple[, 2])
  delta.uniform.range <- c(0, 100)
  grid <- seq(40, 60, by = 4) 
  # grid interval "by = 0.1" is recommended.

  ###  Running the following codes takes more time than CRAN policy
  ###  Please type the following lines without "#" to run the function and to see the results
  #  logprof <- entirelogprofilelikelihood(data = simple, grid = grid, 
  #                                        initial = c(theta.ini, c.ini), data.flux = FALSE, 
  #                                        delta.uniform.range = delta.uniform.range)

  #  plot(grid, logprof, type = "l", 
  #       xlab = expression(bold(Delta)),       
  #       ylab = expression(bold(paste("log L"[prof], "(", Delta, ")"))))
  #  prof <- exp(logprof - max(logprof))  # normalization
  #  plot(grid, prof, type = "l", 
  #       xlab = expression(bold(Delta)),       
  #       ylab = expression(bold(paste("L"[prof], "(", Delta, ")"))))



}


\keyword{methods}
