#data(sysdata, envir=environment())

lineup <- function(m, dataprob=NULL, nulls=NULL) {
  # assume first element is data
  probs <- stats::runif(m)
  n.targets <- length(dataprob)

  if (!is.null(dataprob)) probs[1:n.targets] <- dataprob
  if (!is.null(nulls)) probs[(1+n.targets):m] <- nulls
  sample(m, size=1, prob=1-probs)
}

pickData <- function(m, xp=1, dataprob=NULL, nulls=NULL) {
  probs <- stats::runif(m)
  n.targets = length(dataprob)
  if (!is.null(dataprob)) {
    probs[1:n.targets] <- dataprob
  }
  if (!is.null(nulls)) probs[(n.targets+1):m] <- nulls
#  sample(m, size=1, prob=1-probs)
#  stats::rbinom(1, size=1, prob=f(probs))
  ps <- (1-probs)^xp
  if (all (ps==0)) ps <- rep(1, length(probs))
  stats::rbinom(1, size=1, prob=sum(ps[1:n.targets])/sum(ps))
}

# too deterministic
# function(ps) ps[1] < min(ps[-1]))
# just wrong:
# function(ps) stats::rbinom(1, size=1, prob=(1-ps[1])/(2-ps[1] - min(ps[-1])))

scenario1 <- function(N, K, m = 20, xp=1, target=1) {
# new lineup in each evaluation: new data, new sample of nulls  
  table(replicate(N,  {
    individual <- sum(replicate(K, pickData(m, dataprob=NULL, xp=xp)) %in% target)
    individual
  }))/N
}


scenario2 <- function(N, K, m = 20, xp=1, target=1) {
  # each data evaluated K times, always with different nulls
  table(replicate(N,  {
    n.targets = length(target)
    dataprob <- stats::runif(n.targets)
    individual <- sum(replicate(K, pickData(m, dataprob=dataprob, xp=xp) %in% target))
    individual
  }))/N
}
  

scenario3 <- function(N, K, m=20, xp=1, target=1) {
  # each data evaluated K times, with the same nulls
  table(replicate(N/100,  
    replicate(100, {
      n.targets = length(target)
      dataprob <- stats::runif(n.targets)
      nulls <- stats::runif(m-n.targets)
      
      individual <- sum(replicate(K, pickData(m, dataprob=dataprob, nulls=nulls, xp=xp)) %in% target)
      individual
    })))/N
}

scenario4 <- function(N, K, m=20, xp=1, target=1) {
  # K is vector: length(K) lineups are shown K[i] times to observers
  # all length(K) lineups show the same data, but have different nulls
  
  res <- replicate(N, {
      n.targets = length(target)
      dataprob <- stats::runif(n.targets)
      
      individual <- rep(NA, length(K))
      for (i in 1:length(K)) {
        nulls <- stats::runif(m-n.targets)              
        individual[i] <- sum(replicate(K[i], pickData(m, dataprob=dataprob, nulls=nulls, xp=xp)) %in% target)
      }
      sum(individual)
    })
  table(res)/N
}

#' Bootstrap based p values for visual inference
#' 
#' Simulation based p value to observe x or more picks of the data plot in K evaluations under the assumption that the data plot is consistent with the null hypothesis.
#' We distinguish between three different scenarios:
#' \itemize{
#' \item Scenario I: in each of K evaluations a different data set and a different set of (m-1) null plots is shown.
#' \item Scenario II: in each of K evaluations the same data set but a different set of (m-1) null plots is shown.
#' \item Scenario III: the same lineup, i.e. same data and same set of null plots, is shown to K different observers.
#' }
#' @param x number of observed picks of the data plot
#' @param K number of evaluations
#' @param m size of the lineup
#' @param N MC parameter: number of replicates on which MC probabilities are based. Higher number of replicates will decrease MC variability.
#' @param scenario numeric value, one of 1,2, or 3, indication the type of simulation used: scenario 3 assumes that the same lineup is shown in all K evaluations
#' @param xp exponent used, defaults to 1
#' @param target (vector) of integer values between 1 and m indicating the position(s) of the target plots. Only the number of targets will affect the probabilities.
#' @param upper.tail compute probabilities P(X >= x). Be aware that the use of this parameter is not consistent with the other distribution functions in base. There, a value of P(X > x) is computed for upper.tail=TRUE.
#' @return Vector/data frame. For comparison a p value based on a binomial distribution is provided as well.
#' @export
#' @examples
#' pVsim(15, 20, m=3) # triangle test
pVsim <- function(x, K, m=20, N=10000, scenario=3, xp=1, target=1, upper.tail=TRUE) {
  type <- paste("scenario", scenario, sep="")
  freq <- get(type)(N=N, K=K, m=m, xp=xp, target=target)
  if (upper.tail) {
    sim <- sapply(x, function(y) sum(freq[as.numeric(names(freq)) >= y]))
    return(cbind(x=x, "simulated"=sim, "binom"=1-stats::pbinom(x-1, size=K, prob=1/m)))
  } else {
    sim <- sapply(x, function(y) sum(freq[as.numeric(names(freq)) < y]))
    return(cbind(x=x, "simulated"=sim, "binom"= stats::pbinom(x-1, size=K, prob=1/m)))
  }
}

#' Simulation based density for visual inference
#' 
#' Probablity of observing exactly x picks of the data plot in K evaluations of a lineup of size m.
#' We distinguish between three different scenarios:
#' \itemize{
#' \item Scenario I: in each of K evaluations a different data set and a different set of (m-1) null plots is shown.
#' \item Scenario II: in each of K evaluations the same data set but a different set of (m-1) null plots is shown.
#' \item Scenario III: the same lineup, i.e. same data and same set of null plots, is shown to K different observers.
#' }
#' @param x number of observed picks of the data plot
#' @param K number of evaluations of the same lineup
#' @param m size of the lineup
#' @param N MC parameter: number of replicates on which MC probabilities are based. Higher number of replicates will decrease MC variability.
#' @param scenario numeric value, one of 1, 2, or 3, indicating the type of simulation used: scenario 3 assumes that the same lineup is shown in all K evaluations
#' @param xp exponent used, defaults to 1
#' @param target location of target plot(s). By default 1. If several targets are present, specify vector of target locations.
#' @return simulation based density to observe x picks of the data plot in K evaluation under the assumption that the data plot is consistent with the null hypothesis. For comparison a p value based on a binomial distribution is provided as well.
#' @export
#' @examples
#' dVsim(2, 20, m=3) # triangle test
#' 
#' \dontrun{
#' ## points in red are binomial distribution, black points are for inference
#' ## in lineups using scenario 3
#' require(ggplot2)
#' qplot(x=x, y=scenario3, data=dVsim(0:6,6,m=2)) + 
#'    geom_point(aes(x,y=binom), colour="red") + ylim(c(0,0.5))
#' qplot(x=x, y=scenario3, data=dVsim(0:6,6,m=3)) + 
#'    geom_point(aes(x,y=binom), colour="red") + ylim(c(0,0.5))
#' }
#'    
#' # lineup with two targets: what are the probabilities to identify at least 
#' # one of the targets?
#' dVsim(0:5, K=5, m=20, N=10000, scenario=3, target=1:2)
#' # slight difference between this distribution and the distribution for a 
#' # lineup of size 10 with a single target:
#' dVsim(0:5, K=5, m=10, N=10000, scenario=3, target=1)
dVsim <- function(x, K, m=20, N=10000, scenario=3, xp=1, target=1) {
  argx <- x
  freqs <- data.frame(Var1=0:K)
  for (t in scenario) {
    t <- paste("scenario", t, sep="")
    freq <- data.frame(get(t)(N=N, K=K, m=m, xp=xp, target=target))
    names(freq)[2] <- t
    freqs <- merge(freqs, freq, by="Var1", all=T)
  }
  freq <- freqs
  freq[is.na(freq)] <- 0
  freq$binom <- stats::dbinom(0:K, size=K, prob=length(target)/m)
  names(freq)[1] <- "x"
  subset(freq, x %in% argx)
}

#' Theoretical distribution for lineups under different scenarios
#' 
#' Some more details to be written later
#' @param x vector of the number of picks of the data plot out of K evaluations
#' @param K number of evaluations
#' @param m size of the lineup
#' @param scenario which scenario should be used? 1, 2, or 3?
#' @param type one of "mpfr" or "numeric". Should the result be in arbitrary numeric length or be a numeric? Internally the Rmpfr package is used to get maximal precision.
#' @export
#' @examples 
#' pV(0:5, 5, m=3, scenario=3)
pV <- function(x, K, m, scenario, type="numeric") {
  res <- x
  if (3 %in% scenario) res <- cbind(res, scenario3=pv3(x, K, m, type=type))
  if (1 %in% scenario) res <- cbind(res, scenario1=1-stats::pbinom(x, size=K, prob=1/m)+stats::dbinom(x, size=K, prob=1/m))
  if (2 %in% scenario) res <- cbind(res, scenario2=pv2(x, K, m))

  if (ncol(res) == 2) {
    res <- as.vector(res[,2])
    names(res) <- x
    return(res)
  }
  names(res)[1] <- "x"
  res
}


#' Theoretical density for lineups under different scenarios
#' 
#' Some more details to be written later
#' @param x vector of the number of picks of the data plot out of K evaluations
#' @param K number of evaluations
#' @param m size of the lineup
#' @param scenario which scenario should be used? 1, 2, or 3?
#' @param type one of "mpfr" or "numeric". Should the result be in arbitrary numeric length or be a numeric? Internally the Rmpfr package is used to get maximal precision.
#' @export
#' @examples
#' dV(0:5, 5, m=2, scenario=3)
#' ## compare to 
#' dVsim(0:5, 5, m=2, scenario=3)
#' 
#' require(ggplot2)
#' ## probabilities can be computed without numeric loss for K=50:
#' K <- 50
#' print(qplot(0:K, dV(0:K, K, m=20, scenario=3))); 
#' print(sum(dV(0:K, K, m=20, scenario=3))); 
dV <- function(x, K, m, scenario, type="numeric") {
  res <- x
  if (3 %in% scenario) res <- cbind(res, scenario3=dv3(x, K, m, type=type))
  if (1 %in% scenario) res <- cbind(res, scenario1=stats::dbinom(x, size=K, prob=1/m))
  if (2 %in% scenario) res <- cbind(res, scenario2=dv2(x, K, m))

  if (ncol(res) == 2) {
    res <- as.vector(res[,2])
    names(res) <- x
    return(res)
  }
  names(res)[1] <- "x"
  res
}


#' Theoretical quantile for lineups under different scenarios
#' 
#' Some more details to be written later
#' @param q (vector) of the quantiles for picks of the data plot out of K evaluations
#' @param K number of evaluations
#' @param m size of the lineup
#' @param scenario which scenario should be used? 1, 2, or 3?
#' @param type one of "mpfr" or "numeric". Should the result be in arbitrary numeric length or be a numeric? Internally the Rmpfr package is used to get maximal precision.
#' @export
#' @examples
#' ## get critical values of visual triangle test:
#' qV(q=c(0.95, 0.99), K=c(5,10,15,20, 25, 30), m=3, scenario=1:3)
#' 
#' ## get critical values of full lineup test:
#' qV(q=c(0.95, 0.99), K=c(5,10,15,20, 25, 30), m=20, scenario=3)
qV <- function(q, K, m, scenario, type="numeric") {
  res <- data.frame(expand.grid(q=q, K=K))
  if (1 %in% scenario) {
    # res$scenario1 <- unlist(plyr::alply(res, .margins=1, function(x) stats::qbinom(x$q, size=x$K, prob=1/m)))
    res$scenario1 <- purrr::pmap_dbl(res, function(q, K) stats::qbinom(q, size=K, prob=1/m))
  }
  if (2 %in% scenario) {
    res2 <- qv2(q, K, m)
    names(res2)[3] <- "scenario2"
    res <- merge(res, res2, by=c("q", "K"))
  }
  if (3 %in% scenario) {
    res3 <- qv3(q, K, m)
    names(res3)[3] <- "scenario3"
    res <- merge(res, res3, by=c("q", "K"))
  }
  
  res
}

dv3 <- function(x, K, m, type="numeric") {
  T1 <- function(m) {
    if (m==2) return(expression(1/m*(1/u^2*(log((u+1)/u)*u^2 + u - log(u+1)))))
    if (m==3) return(expression(1/m*(1/u^3*((3*u+1)*log(u) - (u^3-3*u-1)*log((u+1)/u)+ 
                                         (4*u^3-3*u-1)*log((2*u+1)/(2*u)) + 2*u^2 -
                                         3*u*log(2*u) + log(u+1) - log(2*u)))))
  }
  ci <- function(i, K, x) {
    res <- Rmpfr::mpfr(rep(0, length=length(i)), 120)
    idx <- which(i >= K-x)
    res[idx] <- ((-1)^(i[idx]-K+x) * Rmpfr::chooseMpfr(x, i[idx]-K+x))
#    browser()
    res
  }
  Tone <- function(i, m) {
#    load("data/Dis.RData")
#    data(Dis)
    if (i==0) return(1)
    u <- 1
    if (i==1) return(eval(T1(m)))
    (-1)^(i-1)/factorial(i-1)*Dis[[m]][i] #eval(Dks[[i]])
  }
  Tone_saved <- function(i,m) {
    if (i == 0) return(1)
#    data(T1m, package="vinference")
    T1m[[m]][i]
  }
  Tim <- function(i, m) {
    sapply(i, function(j) Tone_saved(j,m))
  }
  Dk <- function(express, name, k) {
    ## returns all k derivatives of express up to k
    res <- list()
    res[[1]] <- express # zeroeth derivative is expression itself
    for (j in 2:(k+1)) res[[j]] <- stats::D(res[[j-1]], name)
    res
  }
  hone <- function (x, K, m) {
#    data(T1m, package="vinference")
    cis <- ci(0:K, K, x)

    #   choose(K, x)*sum(cis*unlist(Tim(0:K, m)))
 #   browser()
    xs <- T1m[[m]][1:K]*Rmpfr::chooseMpfr(K, x)*cis[-1]
    sum(xs) + cis[1]
  }
  
  if (is.null(T1m[[m]])) {
#  if (!(m %in% c(2,3,20))) {
    stop("Not implemented for this lineup size, use bootstrap simulation with dVsim instead.")  
  }
  if (K > length(T1m[[m]])) {
    if (m==3) return(p3(x, K))
    msg <- sprintf("Not implemented for values of K > %d because of large memory needs. Use bootstrap simulation with dVsim instead.", length(T1m[[m]]))
    stop(msg)  
  }
#  Dks <- Dk(T1(m), "u", k=K)
  xs <- lapply(x, hone, K=K, m=m)
  ys <- xs[[1]]
  if (length(xs) > 1)
    for (i in 2:length(xs)) ys <- c(ys, xs[[i]])
  if (type == "mpfr") return(ys)
  else return(as.numeric(ys))
}

qv3 <- function(q, K, m) {
  dframe <- data.frame(expand.grid(q, K))
  names(dframe) <- c("q", "K")
#  require(plyr)
  # res <- plyr::ddply(dframe, c("q", "K"), function(x) {
  #   hs <- cumsum(dv3(x=0:x$K, K=x$K, m=m))
  #   which(hs>=x$q)[1]
  # })
  res <- purrr::pmap_df(dframe, function(q, K) {
    hs <- cumsum(dv3(x=0:K, K=K, m=m))
    data.frame(q = q, K = K, x = which(hs>=q)[1])
  })
  names(res)[3] <- "x"
  res$x[res$x > res$K] <- NA
  res
}
 

dv2 <- function(x,K, m=3) { 
  dv2one <- function(x, K, m=m) {
    g <- function(q) {
      res <- q*((q+2)*log(q+2)-2*(q+1)*log(q+1)+q*log(q))
      res[q==0] <- 0
      res
    }
    if (!(m %in% c(3,20))) stop("Function not implemented for lineup sizes <> 3. Use bootstrap based density estimation in dVsim instead.")
    
    f <- function(q, K, x) choose(K,x)*g(q)^x*(1-g(q))^(K-x)
    stats::integrate(f, 0,1,K=K,x=x)$value
  }
  unlist(sapply(x, dv2one, K=K, m=m))
}

pv2 <- function(x,K, m=3) { 
  hdone <- function(x1, K, m) {
    sum(dv2(x1:K, K, m))
  }
  sapply(x, hdone, K=K, m=m)
}

qv2 <- function (q, K, m=3) {
  dframe <- data.frame(expand.grid(q, K))
  names(dframe) <- c("q", "K")
#  require(plyr)
  # res <- plyr::ddply(dframe, c("q", "K"), function(x) {
  #   hs <- cumsum(dv2(x=0:x$K, K=x$K, m=m))
  #   which(hs>=x$q)[1]
  # })
  res <- purrr::pmap_df(dframe, function(q, K) {
    hs <- cumsum(dv2(x=0:K, K=K, m=m))
    data.frame(q = q, K = K, x = which(hs>=q)[1])
  })
  names(res)[3] <- "x"
  res$x[res$x > res$K] <- NA
  res
}



pv3 <- function(x, K, m, type="numeric") {
  hdone <- function(x1, K, m) {
    sum(dv3(x1:K, K, m))
  }
  sapply(x, hdone, K=K, m=m)
}

