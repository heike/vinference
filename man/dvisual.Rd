\name{dvisual}
\alias{dvisual}
\title{Bootstrap based density for visual inference}
\usage{
  dvisual(x, K, m = 20, N = 10000, type = "scenario4")
}
\arguments{
  \item{x}{number of observed picks of the data plot}

  \item{K}{number of evaluations of the same lineup}

  \item{m}{size of the lineup}

  \item{N}{MC parameter: number of replicates on which MC
  probabilities are based. Higher number of replicates will
  decrease MC variability.}

  \item{type}{type of simulation used: scenario 4 assumes
  that the same lineup is shown in all K evaluations}
}
\value{
  simulation based density to observe x picks of the data
  plot in K evaluation under the assumption that the data
  plot is consistent with the null hypothesis. For
  comparison a p value based on a binomial distribution is
  provided as well.
}
\description{
  Probablity of observing exactly x picks of the data plot
  in K evaluations of the same lineup of size m.
}
\examples{
dvisual(2, 20, m=3) # triangle test
}

