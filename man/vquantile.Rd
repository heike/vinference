\name{vquantile}
\alias{vquantile}
\title{Quantile of the analytically derived quantiles for visual inference}
\usage{
  vquantile(q, K, m,
    type = c("scenario1", "scenario2", "scenario3"))
}
\arguments{
  \item{q}{(vector) of quantiles}

  \item{K}{number of evaluations}

  \item{m}{lineup size, currently only m=2 and 3 are
  treated analytically. Use simulation within dvisual to
  get to other values for m}

  \item{type}{which scenario was used? One of scenario1,
  scenario2, scenario3}
}
\description{
  Quantile of the analytically derived quantiles for visual
  inference
}
\examples{
## get critical values of visual triangle test:
hquantile(q=c(0.95, 0.99), K=c(5,10,15,20, 25, 30), m=3)

## get critical values of full lineup test:
vquantile(q=c(0.95, 0.99), K=c(5,10,15,20, 25, 30), m=20)
}

