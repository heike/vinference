\name{dvisual}
\alias{dvisual}
\title{Bootstrap based density for visual inference}
\usage{
  dvisual(x, K, m = 20, N = 10000, type = "scenario3",
    xp = 1)
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

  \item{xp}{exponent used, defaults to 1}
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
  in K evaluations of a lineup of size m. We distinguish
  between three different scenarios: \itemize{ \item
  Scenario I: in each of K evaluations a different data set
  and a different set of (m-1) null plots is shown. \item
  Scenario II: in each of K evaluations the same data set
  but a different set of (m-1) null plots is shown. \item
  Scenario III: the same lineup, i.e. same data and same
  set of null plots, is shown to K different observers. }
}
\examples{
dvisual(2, 20, m=3) # triangle test
qplot(x=x, y=scenario3, data=dvisual(0:6,6,m=2)) + geom_point(aes(x,y=binom), colour="red") + ylim(c(0,0.5))
qplot(x=x, y=scenario3, data=dvisual(0:6,6,m=3)) + geom_point(aes(x,y=binom), colour="red") + ylim(c(0,0.5))
}

