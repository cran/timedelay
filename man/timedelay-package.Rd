\name{timedelay}
\alias{timedelay-package}
\alias{timedelay}
\docType{package}
\title{Time Delay Estimation for Stochastic Time Series of Gravitationally Lensed Quasars}
\description{
The R package \pkg{timedelay} provides a toolbox to estimate the time delay between the brightness time series of gravitationally lensed quasar images via Bayesian and profile likelihood approaches. The model is based on a state-space representation for  irregularly observed time series data generated from a latent continuous-time Ornstein-Uhlenbeck process. Our Bayesian method adopts scientifically motivated hyper-prior distributions and a Metropoli-Hastings within Gibbs sampler, producing posterior samples of the model parameters that include the time delay. A profile likelihood of the time delay is a simple approximation to the marginal posterior distribution of the time delay. Both Bayesian and profile likelihood approaches complement each other, producing almost identical results; the Bayesian way is more principled but the profile likelihood is easier to be implemented.
}
\details{
  \tabular{ll}{
    Package: \tab timedelay\cr
    Type: \tab Package\cr
    Version: \tab 1.0.0\cr
    Date: \tab 2015-05-25\cr
    License: \tab GPL-2\cr
    Main functions: \tab \code{\link{bayesian}}, \code{\link{logprofilelikelihood}}, \code{\link{entirelogprofilelikelihood}}\cr
  }

}
\author{
  Hyungsuk Tak, Kaisey Mandel, David A. van Dyk, Vinay L. Kashyap, Xiao-Li Meng, and Aneta Siemiginowska

  Maintainer: Hyungsuk Tak <hyungsuk.tak@gmail.com>
}
\references{
  Hyungsuk Tak, Kaisey Mandel, David A. van Dyk, Vinay L. Kashyap, Xiao-Li Meng, and Aneta Siemiginowska (in progress). Bayesian and Profile Likelihood Approaches to Time Delay Estimation for Stochastic Time Series of Gravitationally Lensed Quasars}

\examples{

  # Loading datasets
  data(simple)
  head(simple)

  ###############################################
  # Time delay estimation via Bayesian approach #
  ###############################################

  output <- bayesian(dat = simple, data.flux = FALSE, theta.ini = c(0, 0.03, 100), 
                     delta.ini = 50, delta.uniform.range = c(0, 100), delta.proposal.scale = 0.3, 
                     tau.proposal.scale = 2, tau.prior.shape = 1, tau.prior.scale = 1,
                     sigma.prior.shape = 1, sigma.prior.scale = 2 * 10^(-7), asis = TRUE,
                     sample.size = 50, warmingup.size = 50)

  names(output)
  hist(output$delta, 30, xlab = expression(bold(Delta)), 
       main = expression(bold(paste("Histogram of posterior ", Delta))))
  plot(output$delta, ylab = expression(bold(Delta)), xlab = "Iteration",
       main = expression(bold(paste("Traceplot of posterior ", Delta))), type = "l")
  acf(output$delta, main = expression(bold(paste("ACF of posterior ", Delta))))

  ### Argument description

  # data.flux: "True" if data are recorded on flux scale or "FALSE" if data are on magnitude scale.
  # theta.ini: Initial values of theta = (mu, sigma, tau) for MCMC.
  # delta.ini: Initial values of the time delay for MCMC.
  # delta.uniform.range: The range of the Uniform prior distribution for the time delay.
  #                      The feasible entire support is 
  #                      c(min(simple[, 1]) - max(simple[, 1]), max(simple[, 1]) - min(simple[, 1]))

  # delta.proposal.scale: The proposal scale of the Metropolis step for the time delay.
  # tau.proposal.scale: The proposal scale of the Metropolis-Hastings step for tau.

  # tau.prior.shape: The shape parameter of the Inverse-Gamma hyper-prior distribution for tau. 
  # tau.prior.scale: The scale parameter of the Inverse-Gamma hyper-prior distribution for tau. 
  # sigma.prior.shape: The shape parameter of 
  #                    the Inverse-Gamma hyper-prior distribution for sigma^2.
  # sigma.prior.scale: The scale parameter of 
  #                    the Inverse-Gamma hyper-prior distribution for sigma^2.

  # asis: "TRUE" if we use the ancillarity-sufficiency interweaving strategy (ASIS) 
  #       for c (always recommended)

  # adaptive.freqeuncy: The adaptive MCMC is applied for every specified frequency. 
  #                     If it is specified as 500, 
  #                     the adaptive MCMC is applied to every 500th iterstion.
  # adaptive.delta: "TRUE" if we use the adaptive MCMC for the time delay.
  # adaptive.delta.factor: The factor, exp(adaptive.delta.factor) or exp(-adaptive.delta.factor), 
  #                        multiplied to the proposal scale of the time delay for adaptive MCMC.
  # adaptive.delta.acceptance.rate.upper.bound: We multiply exp(adaptive.delta.factor) 
  #                                             to the proposal scale of the time delay 
  #                                             if the acceptance rate of the time delay 
  #                                             is greater than this number
  # adaptive.delta.acceptance.rate.lower.bound: We multiply exp(-adaptive.delta.factor) 
  #                                             to the proposal scale of the time delay 
  #                                             if the acceptance rate of the time delay 
  #                                             is smaller than this number

  # adaptive.tau: "TRUE" if we use the adaptive MCMC for tau.
  # adaptive.tau.factor: The factor, exp(adaptive.tau.factor) or exp(-adaptive.tau.factor), 
  #                      multiplied to the proposal scale of tau for adaptive MCMC.
  # adaptive.tau.acceptance.rate.upper.bound: We multiply exp(adaptive.tau.factor) 
  #                                           to the proposal scale of tau
  #                                           if the acceptance rate of tau
  #                                           is greater than this number
  # adaptive.tau.acceptance.rate.lower.bound: We multiply exp(-adaptive.tau.factor) 
  #                                           to the proposal scale of tau
  #                                           if the acceptance rate of tau
  #                                           is smaller than this number

  # tempered.transition: "TRUE" if we use the tempered transition.
  #                      We recommend not using the adaptive MCMC for the time delay 
  #                      together with the tempered transition,
  #                      i.e., "adaptive.delta = FALSE" if "tempered.transition = TRUE"
  # number.rungs: The number of rungs in the temperature ladder
  # temperature.base: The temerature of the i-th rung is (temperature.base)^i

  # sample.size: The number of posterior samples for each parameter
  # warmingup.size: The number of burn-ins

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

  ### Argument description

  # Delta: The log profile likelihood value at this number.
  # initial: The initial values of the other model parameters (mu, sigma, tau, c)
  # data: The data set   
  # data.flux: "True" if data are recorded on flux scale or "FALSE" if data are on magnitude scale.
  # delta.uniform.range: The range of the Uniform prior distribution for the time delay.

  # This function "logprofilelikelihood" can be used to find the maximum profile likelihood estimate
  # and the inverse of the negative Hessian at the mode.
  # See the example in the manual by typing "?logprofilelikelihood" at the R console.

  ###### The entire profile likelihood values on the grid of values of the time delay.

  grid <- seq(45, 55, by = 5)
  # grid interval "by = 0.1" is recommended for accuracy.

  ###  Running the following codes takes more time than CRAN policy
  ###  Please type the following lines without "#" to run the function and to see the results
  #  logprof <- entirelogprofilelikelihood(data = simple, grid = grid, 
  #                                        initial = c(theta.ini, c.ini), data.flux = FALSE, 
  #                                        delta.uniform.range = delta.uniform.range)
  #  plot(grid, logprof, type = "l", 
  #       xlab = expression(bold(Delta)),       
  #       ylab = expression(bold(paste("log L"[prof], "(", Delta, ")"))))
  #  prof <- exp(logprof - max(logprof))
  #  plot(grid, prof, type = "l", 
  #       xlab = expression(bold(Delta)),       
  #       ylab = expression(bold(paste("L"[prof], "(", Delta, ")"))))

  ### Argument description

  # data: The data set   
  # grid: the vector of grid values of the time delay 
  #       on which the log profile likelihood values are calculated.
  # initial: The initial values of the other model parameters (mu, sigma, tau, c)
  # data.flux: "True" if data are recorded on flux scale or "FALSE" if data are on magnitude scale.
  # delta.uniform.range: The range of the Uniform prior distribution for the time delay.

}

\keyword{package}