#' Estimate alpha of a symmetric Dirichlet distribution
#' 
#' The estimation of the rate parameter alpha is done with Maximum Likelihood Estimation.
#' @param p (matrix) of probabilities. Number of columns is interpreted as dimension of the Dirichlet distribution
#' Number of rows is repetitions. 
#' @param weight vector
#' @param eps error allowed in finding root
#' @return rate parameter alpha of symmetric (flat) Dirichlet distribution
#' @importFrom stats rgamma uniroot
#' @export
#' @examples
#' x <- matrix(rgamma(n=2000, shape=.1), ncol=20)
#' plot(density(x[,1]))
#' x <- x/rowSums(x)
#' alpha_ml(x)
#' shapes <- seq(0.01, 0.75, by = 0.01)
#' ests <- sapply(shapes, function(s) {
#'   x <- matrix(rgamma(n=2000, shape=s), ncol=20)
#'   alpha_ml(x, eps=10^-20)
#' })
#' plot(shapes,ests)
#' abline(a=0,b=1)
alpha_ml <- function(p, weight=NULL, eps=10^(-7)) {
  # matrix p
  if (is.null(dim(p))) {
    m <- length(p)
    n <- 1
    p <- matrix(p, nrow=1)
  } else {
    n <- dim(p)[1]
    m <- dim(p)[2]
  }
  if (is.null(weight)) weight <- rep(1, n)
  else n <- sum(weight)
  
  weight <- weight/sum(weight)
  ps <- p
  ps[p < eps] <- eps # make sure we don't take the log of 0
  ps <- ps/rowSums(ps)
  logp <- sum(weight*rowSums(log(ps)))/m
  ml <- function(alpha) {
    digamma(alpha) - digamma(alpha*m) - logp
  }
  # find alpha such that digamma(alpha) - digamma(alpha*m) = logp
  alpha <- uniroot(f=ml, interval=c(eps,10))
  
  
  alpha$root
}
