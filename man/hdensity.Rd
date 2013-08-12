\name{hdensity}
\alias{hdensity}
\title{Theoretical density for lineups under scenario 3 for m = 2,3}
\usage{
  hdensity(x, K, m, type = "numeric")
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
hdensity(0:5, 5, m=2)
## compare to
dvisual(0:5, 5, m=2)

## test how many K can be computed properly
for (K in 10:50) {
  print(K);
  print(qplot(0:K, hdensity(0:K, K, m=20)));
  print(sum(hdensity(0:K, K, m=20)));
  scan()
}
}

