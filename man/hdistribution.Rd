\name{hdistribution}
\alias{hdistribution}
\title{Theoretical distribution for lineups under scenario 3 for m = 2, 3, and 20}
\usage{
  hdistribution(x, K, m, type = "numeric")
}
\arguments{
  \item{x}{vector of the number of picks of the data plot
  out of K evaluations}

  \item{K}{number of evaluations}

  \item{m}{size of the lineup}

  \item{type}{one of "mpfr" or "numeric". Should the result
  be in arbitrary numeric length or be a numeric?
  Internally the Rmpfr package is used to get maximal
  precision.}
}
\description{
  Some more details to be written later
}
\examples{
hdistribution(0:5, 5, m=3)
}

