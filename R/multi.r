#' @include pvalues.r
dvismulti1 <- function(K, k, m=20, N=5000) {
  stopifnot((K==length(k)) | (length(k)==1))
  k <- rep(k, length=K)
  res <- replicate(N, {
    # everybody is shown a different lineup
 #   require(plyr)
    success <- purrr::map_dbl(1:K, function(i) {
      lp <- stats::runif(m)
      sum(sample(1:m, size=k[i], replace=FALSE, prob=lp)==1)
    })
    # success <- plyr::ldply(1:K, )
    sum(success)
  })
  res <- factor(res, levels=0:K)
  table(res)/N
}

dvismulti2 <- function(K, k, m=20, N=5000) {
  stopifnot((K==length(k)) | (length(k)==1))
  k <- rep(k, length=K)
  res <- replicate(N, {
    # same data, different nulls
 #   require(plyr)
    first <- stats::runif(1)
    success <- purrr::map_dbl(1:K, function(i) {
      lp <- c(first, stats::runif(m-1))
      sum(sample(1:m, size=k[i], replace=FALSE, prob=lp)==1)
    })
    # success <- plyr::ldply(1:K, function(i) {
    #   lp <- c(first, runif(m-1))
    #   sum(sample(1:m, size=k[i], replace=FALSE, prob=lp)==1)
    # })
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
    lp <- stats::runif(m)
 #   require(plyr)
    success <- purrr::map_dbl(1:K, function(i) {
      sum(sample(1:m, size=k[i], replace=FALSE, prob=lp)==1)
    })
    # success <- plyr::ldply(1:K, function(i) {
    #   sum(sample(1:m, size=k[i], replace=FALSE, prob=lp)==1)
    # })
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
#' @name vmulti
#' @aliases dmulti pmulti qmulti 
#' @param x (vector) of the number of observed data picks
#' @param q (vector) of the number of observed data picks
#' @param K integer value, specifying the number of participants/evaluations 
#' @param k (vector) of integer values of the number of observed picks in a multiple choice lineup evaluation
#' @param m lineup size, defaults to m=20
#' @param type character, one of "scenario1", "scenario2", or "scenario3"  
#' @param N integer value, number of simulations. 
#' @return list consisting of two named items: 
#' \itemize{
#' \item dmulti: density estimates and vector k 
#' \item pmulti: distribution estimates and vector k 
#' \item qmulti: quantile estimates and vector k 
#' }
#' @examples
#' k=rpois(5,lambda=1)+1
#' m=20
#' dmulti(0:5,K=5,k=k, type="scenario1", m) 
#' ## compare to Poisson Binomial:
#' require(poibin)
#' dpoibin(0:5, pp=k/m)
#' 
#' qmulti(c(0.95, 0.99), K=5, k=2)
#' 
#' \dontrun{
#' (k <- rpois(5, lambda=1.5)+1)
#' reps1 <- plyr::ldply(1:10, function(x) dmulti(0:5, 5, k, type="scenario1")$density)
#' reps2 <- plyr::ldply(1:10, function(x) dmulti(0:5, 5, k, type="scenario2")$density)
#' reps3 <- plyr::ldply(1:10, function(x) dmulti(0:5, 5, k, type="scenario3")$density)
#' reps1$type <- "I"
#' reps2$type <- "II"
#' reps3$type <- "III"
#' reps <- rbind(reps1, reps2, reps3)
#' library(reshape2)
#' library(RColorBrewer)
#' cols <- brewer.pal(8, "Paired")
#' 
#' mr <- melt(reps, measure.vars=1:6)
#' mr$variable <- as.numeric(gsub("V", "", mr$variable))
#' mr$type <- factor(mr$type)
#' mr$offset <- with(mr, 0.1*c(0,-1,1)[as.numeric(type)])
#' qplot(variable+offset, value, data=mr, alpha=I(0.75), colour=type) + 
#'   geom_point(aes(x=1:6, y=dpoibin(0:5, pp=k*1/m)), pch=1, size=4, 
#'     inherit.aes=FALSE)+
#'   xlab("x")+ylab("P(X=x)")+ggtitle(paste(k, collapse=",")) + 
#'   theme_bw() + scale_colour_manual(values=cols[-c(1,3,5,6,7)])
#'  }


NULL

#' @rdname vmulti
#' @aliases dmulti pmulti qmulti 
#' @export
dmulti<- function(x, K, k, m=20, type="scenario3", N=5000) {
  density <- switch(type, scenario1 = dvismulti1(K,k,m,N),
         scenario2 = dvismulti2(K,k,m,N),
         scenario3 = dvismulti3(K,k,m,N))
  res <-list(density=as.vector(density)[x+1], k=k)
  
  res
}

#' @rdname vmulti
#' @aliases dmulti pmulti qmulti 
#' @export
pmulti<- function(x, K, k, m=20, type="scenario3", N=5000) {
  density <- switch(type, scenario1 = dvismulti1(K,k,m,N),
                    scenario2 = dvismulti2(K,k,m,N),
                    scenario3 = dvismulti3(K,k,m,N))
  res <-list(distribution=as.vector(cumsum(density))[x+1], k=k)
  
  res
}

#' @rdname vmulti
#' @aliases dmulti pmulti qmulti 
#' @export
qmulti<- function(q, K, k, m=20, type="scenario3", N=5000) {
  density <- switch(type, scenario1 = dvismulti1(K,k,m,N),
                    scenario2 = dvismulti2(K,k,m,N),
                    scenario3 = dvismulti3(K,k,m,N))

  res <-list(quantile=cumsum(as.vector(density)), k=k)
  # res$quantile <- plyr::laply(q, function(qq) min(which(res$quantile >= qq)-1))
  res$quantile <- purrr::map_dbl(q, function(qq) min(which(res$quantile >= qq)-1))
  names(res$quantile) <- q
  
  res
}