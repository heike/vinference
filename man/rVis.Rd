% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_alpha.R
\name{rVis}
\alias{rVis}
\title{Functions to produce visual estimate of alpha
Simulate results from a lineup evaluation experiment using the Dirichlet-Multinomial model}
\usage{
rVis(N = 50, K = 22, m = 20, alpha, scenario = 3)
}
\arguments{
\item{N}{Number of lineups to simulate}

\item{K}{The total number of null panel selections (or, in a Rorschach lineup, the total number of evaluations)}

\item{m}{The number of panels in the lineup}

\item{alpha}{The (scalar) symmetric Dirichlet parameter which is related to the number of interesting panels}

\item{scenario}{which lineup administration is used? scenario 1 and 3 are implemented}
}
\value{
Matrix of dimension m by N
}
\description{
This function returns panel selection counts simulated from the Dirichlet-
Multinomial model; that is, the result is a \eqn{m \times N} matrix of panel selection counts.
}
\examples{
rVis(alpha = .5, m = 20, K = 30, N = 5, scenario = 3)
rVis(alpha = .5, m = 20, K = 30, N = 5, scenario = 1)
}
