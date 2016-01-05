\name{simple}

\docType{data}
\alias{simple}

\title{
  Simulated simple data of a doubly-lensed quasar
}

\description{
  We simulated a simple data set of a doubly-lensed quasar with generative parameter values (Delta, beta0, mu, sigma, tau) equal to (50, 2, 0, 0.03, 100). The number of epochs is 78 and measurement errors are the same to 0.005.
}

\usage{data(simple)}

\format{
  A simple data set of a doubly-lensed quasar:
  \describe{
    \item{\code{time}}{observation times}
    \item{\code{x}}{magnitudes of light A}
    \item{\code{se.x}}{measurement errors of x}
    \item{\code{y}}{magnitudes of light B}
    \item{\code{se.y}}{measurement errors of y}
  }
}

\examples{
  data(simple)
}

\keyword{datasets}
