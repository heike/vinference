\name{pvisual}
\alias{pvisual}
\title{Bootstrap based p values for visual inference}
\usage{
  pvisual(x, K, m = 20, N = 10000, type = "scenario3",
    upper.tail = TRUE)
}
\arguments{
  \item{x}{number of observed picks of the data plot}

  \item{K}{number of evaluations of the same lineup}

  \item{m}{size of the lineup}

  \item{N}{MC parameter: number of replicates on which MC
  probabilities are based. Higher number of replicates will
  decrease MC variability.}

  \item{type}{type of simulation used: scenario 3 assumes
  that the same lineup is shown in all K evaluations}

  \item{upper.tail}{compute probabilities P(X >= x). Be
  aware that the use of this parameter is not consistent
  with the other distribution functions in base. There, a
  value of P(X > x) is computed for upper.tail=TRUE.}
}
\value{
  simulation based p value to observe x or more picks of
  the data plot in K evaluation under the assumption that
  the data plot is consistent with the null hypothesis. For
  comparison a p value based on a binomial distribution is
  provided as well.
}
\description{
  p value of observing at least x picks of the data plot in
  K evaluations of the same lineup of size m.
}
\examples{
pvisual(15, 20, m=3) # triangle test
}

