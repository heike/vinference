\name{hquantile}
\alias{hquantile}
\title{Quantile of the analytically derived density for visual inference according to scenario 3}
\usage{
  hquantile(q, K, m)
}
\arguments{
  \item{q}{(vector) of quantiles}

  \item{K}{number of evaluations}

  \item{m}{lineup size, currently only m=2 and 3 are
  treated analytically. Use simulation within dvisual to
  get to other values for m}
}
\description{
  Quantile of the analytically derived density for visual
  inference according to scenario 3
}
\examples{
## get critical values of visual triangle test:
hquantile(q=c(0.95, 0.99), K=c(5,10,15,20), m=3)
}

