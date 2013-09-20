dvismult1 <- function(K, k, m=20, N=5000) {
  stopifnot((K==length(k)) | (length(k)==1))
  k <- rep(k, length=K)
  res <- replicate(N, {
    # everybody is shown a different lineup
    require(plyr)
    success <- ldply(1:K, function(i) {
      lp <- runif(m)
      sum(sample(1:20, size=k[i], replace=FALSE, prob=lp)==1)
    })
    sum(success)
  })
  res <- factor(res, levels=0:K)
  table(res)/N
}

dvismult2 <- function(K, k, m=20, N=5000) {
  stopifnot((K==length(k)) | (length(k)==1))
  k <- rep(k, length=K)
  res <- replicate(N, {
    # same data, different nulls
    require(plyr)
    first <- runif(1)
    success <- ldply(1:K, function(i) {
      lp <- c(first, runif(m-1))
      sum(sample(1:20, size=k[i], replace=FALSE, prob=lp)==1)
    })
    sum(success)
  })
  res <- factor(res, levels=0:K)
  table(res)/N
}
dvismulti3 <- function(K, k, m=20, N=5000) {
  stopifnot((K==length(k)) | (length(k)==1))
  k <- rep(k, length=K)
  res <- replicate(N, {
    # everybody is shown the same lineup
    lp <- runif(m)
    require(plyr)
    success <- ldply(1:K, function(i) {
      sum(sample(1:20, size=k[i], replace=FALSE, prob=lp)==1)
    })
    sum(success)
  })
  res <- factor(res, levels=0:K)
  table(res)/N
}

#' Simulation based density and distribution functions of visual inference under a multiple choice lineup
#' 
#' Simulation based p value to observe x or more picks of the data plot in K evaluations of a multiple choice lineup test  under the assumption that the data plot is consistent with the null hypothesis.
#' We distinguish between three different scenarios:
#' \itemize{
#' \item Scenario I: in each of K evaluations a different data set and a different set of (m-1) null plots is shown.
#' \item Scenario II: in each of K evaluations the same data set but a different set of (m-1) null plots is shown.
#' \item Scenario III: the same lineup, i.e. same data and same set of null plots, is shown to K different observers.
#' }
#' @param x (vector) of the number of observed data picks
#' @param K integer value, specifying the number of participants/evaluations 
#' @param k (vector) of integer values of the number of observed picks in a multiple choice lineup evaluation
#' @param m lineup size, defaults to m=20
#' @param type character, one of "scenario1", "scenario2", or "scenario3"  
#' @param N integer value, number of simulations. 
#' @return list of two named items: density estimates and vector k 
#' @export
#' @examples
#' dmulti(0:5,5, k=rpois(5,lambda=1)+1)
dmulti<- function(x, K, k, m=20, type="scenario3", N=5000) {
  density <- switch(type, scenario1 = dvismulti1(K,k,m,N),
         scenario2 = dvismulti2(K,k,m,N),
         scenario3 = dvismulti3(K,k,m,N))
  res <-list(density=as.vector(density)[x+1], k=k)
  
  res
}