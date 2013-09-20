\name{dmulti}
\alias{dmulti}
\title{Simulation based density and distribution functions of visual inference under a multiple choice lineup}
\usage{
  dmulti(x, K, k, m = 20, type = "scenario3", N = 5000)
}
\arguments{
  \item{x}{(vector) of the number of observed data picks}

  \item{K}{integer value, specifying the number of
  participants/evaluations}

  \item{k}{(vector) of integer values of the number of
  observed picks in a multiple choice lineup evaluation}

  \item{m}{lineup size, defaults to m=20}

  \item{type}{character, one of "scenario1", "scenario2",
  or "scenario3"}

  \item{N}{integer value, number of simulations.}
}
\value{
  list of two named items: density estimates and vector k
}
\description{
  Simulation based p value to observe x or more picks of
  the data plot in K evaluations of a multiple choice
  lineup test under the assumption that the data plot is
  consistent with the null hypothesis. We distinguish
  between three different scenarios: \itemize{ \item
  Scenario I: in each of K evaluations a different data set
  and a different set of (m-1) null plots is shown. \item
  Scenario II: in each of K evaluations the same data set
  but a different set of (m-1) null plots is shown. \item
  Scenario III: the same lineup, i.e. same data and same
  set of null plots, is shown to K different observers. }
}
\examples{
dmulti(0:5,5, k=rpois(5,lambda=1)+1)
}

