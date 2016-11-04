\name{entirelogprofilelikelihood}			
\alias{entirelogprofilelikelihood}
\title{Calculating the entire profilel likelihood curve over the given grid values of the time delay}
\description{
  \code{entirelogprofilelikelihood} calculates the entire profilel likelihood curve over the given grid values of the time delay. 
}

\usage{entirelogprofilelikelihood(data_lcA, data_lcB, grid, 
                                  initial, data.flux, 
                                  delta.uniform.range, micro)}

\arguments{
  \item{data_lcA}{
    A (\emph{n} by 3) matrix; the first column has n observation times of light curve A, the second column has n flux (or magnitude) values of light curve A, the third column has n measurement errors of light curve A.
  }
  \item{data_lcB}{
    A (\emph{w} by 3) matrix; the first column has w observation times of light curve B, the second column has w flux (or magnitude) values of light curve B, the third column has w measurement errors of light curve B.
  }
  \item{grid}{
    A vector containing values of the time delay on which the profile likelihood values are calculated. We recommend using the grid interval equal to 0.1.
  }
  \item{initial}{
    The initial values of the other model parameters (mu, log(sigma), log(tau), beta). We take log on sigma and tau for numerical stability.
  }
  \item{data.flux}{
    "True" if data are recorded on flux scale or "FALSE" if data are on magnitude scale.
  }
  \item{delta.uniform.range}{
    The range of the Uniform prior distribution for the time delay. The feasible entire support is c(min(simple[, 1]) - max(simple[, 1]), max(simple[, 1]) - min(simple[, 1])).
  }
  \item{micro}{
    It determines the order of a polynomial regression model that accounts for the difference between microlensing trends. Default is 3. When zero is assigned, the Bayesian model fits a curve-shifted model.
  }

}

\details{
  The function \code{entirelogprofilelikelihood} is used to obtain the entire profile likelihood curve over the given grid values of the time delay. 
}

\value{
  The outcome of \code{entirelogprofilelikelihood} is the values of the log profile likelihood function over the given grid values of the time delay.
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
  head(simple)

  ################################################
  # Time delay estimation via profile likelihood #
  ################################################


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

  ###### The entire profile likelihood values on the grid of values of the time delay.

  # Cubic microlensing model
  ti1 <- lcB[, 1]
  ti2 <- lcB[, 1]^2
  ti3 <- lcB[, 1]^3
  ss <- lm(lcB[, 2] - mean(lcA[, 2]) ~ ti1 + ti2 + ti3)

  initial <- c(mean(lcA[, 2]), log(0.01), log(200), ss$coefficients)
  delta.uniform.range <- c(0, 100)
  grid <- seq(0, 100, by = 0.1) 
  # grid interval "by = 0.1" is recommended,
  # but users can set a finer grid of values of the time delay.

  ###  Running the following codes takes more time than CRAN policy
  ###  Please type the following lines without "#" to run the function and to see the results

  #  logprof <- entirelogprofilelikelihood(data_lcA = lcA, data_lcB = lcB, grid = grid, 
  #                                        initial = initial, data.flux = FALSE, 
  #                                        delta.uniform.range = delta.uniform.range, micro = 3)

  #  plot(grid, logprof, type = "l", 
  #       xlab = expression(bold(Delta)),       
  #       ylab = expression(bold(paste("log L"[prof], "(", Delta, ")"))))
  #  prof <- exp(logprof - max(logprof))  # normalization
  #  plot(grid, prof, type = "l", 
  #       xlab = expression(bold(Delta)),       
  #       ylab = expression(bold(paste("L"[prof], "(", Delta, ")"))))



}


\keyword{methods}
