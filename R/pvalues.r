lineup <- function(m, dataprob=NULL, nulls=NULL) {
  # assume first element is data
  probs <- runif(m)
  if (!is.null(dataprob)) probs[1] <- dataprob
  if (!is.null(nulls)) probs[2:(length(nulls)+1)] <- nulls
  sample(m, size=1, prob=1-probs)
}

scenario1 <- function(N, K, m = 20) {
# new lineup in each evaluation: new data, new sample of nulls  
  table(replicate(N,  {
    individual <- sum(replicate(K, lineup(m, dataprob=NULL)) ==1)
    individual
  }))/N
}


scenario2 <- function(N, K, m = 20) {
  # each data evaluated K times, always with different nulls
  table(replicate(N,  {
    dataprob <- runif(1)
    individual <- sum(replicate(K, lineup(m, dataprob=dataprob))==1)
    individual
  }))/N
}
  

scenario3 <- function(N, K, m=20) {
  # each data evaluated K times, with the same nulls
  table(replicate(N/100,  
    replicate(100, {
      dataprob <- runif(1)
      nulls <- runif(m-1)
      
      individual <- sum(replicate(K, lineup(m, dataprob=dataprob, nulls=nulls))==1)
      individual
    })))/N
}

#' Bootstrap based p values for visual inference
#' 
#' p value of observing at least x picks of the data plot in K evaluations of the same lineup of size m.
#' @param x number of observed picks of the data plot
#' @param K number of evaluations of the same lineup
#' @param m size of the lineup
#' @param N MC parameter: number of replicates on which MC probabilities are based. Higher number of replicates will decrease MC variability.
#' @param type type of simulation used: scenario 3 assumes that the same lineup is shown in all K evaluations
#' @param upper.tail compute probabilities P(X >= x). Be aware that the use of this parameter is not consistent with the other distribution functions in base. There, a value of P(X > x) is computed for upper.tail=TRUE.
#' @return simulation based p value to observe x or more picks of the data plot in K evaluation under the assumption that the data plot is consistent with the null hypothesis. For comparison a p value based on a binomial distribution is provided as well.
#' @export
#' @examples
#' pvisual(15, 20, m=3) # triangle test
pvisual <- function(x, K, m=20, N=10000, type="scenario3", upper.tail=TRUE) {
  freq <- get(type)(N=N, K=K, m=m)
  if (upper.tail) {
    sim <- sapply(x, function(y) sum(freq[as.numeric(names(freq)) >= y]))
    return(cbind(x=x, "simulated"=sim, "binom"=1-pbinom(x-1, size=K, prob=1/m)))
  } else {
    sim <- sapply(x, function(y) sum(freq[as.numeric(names(freq)) < y]))
    return(cbind(x=x, "simulated"=sim, "binom"= pbinom(x-1, size=K, prob=1/m)))
  }
}

#' Bootstrap based density for visual inference
#' 
#' Probablity of observing exactly x picks of the data plot in K evaluations of the same lineup of size m.
#' @param x number of observed picks of the data plot
#' @param K number of evaluations of the same lineup
#' @param m size of the lineup
#' @param N MC parameter: number of replicates on which MC probabilities are based. Higher number of replicates will decrease MC variability.
#' @param type type of simulation used: scenario 3 assumes that the same lineup is shown in all K evaluations
#' @return simulation based density to observe x picks of the data plot in K evaluation under the assumption that the data plot is consistent with the null hypothesis. For comparison a p value based on a binomial distribution is provided as well.
#' @export
#' @examples
#' dvisual(2, 20, m=3) # triangle test
#' qplot(x=x, y=simulated, data=data.frame(dvisual(0:6,6,m=2))) + geom_point(aes(x,y=binom), colour="red") + ylim(c(0,0.5))
#' qplot(x=x, y=simulated, data=data.frame(dvisual(0:6,6,m=3))) + geom_point(aes(x,y=binom), colour="red") + ylim(c(0,0.5))
dvisual <- function(x, K, m=20, N=10000, type="scenario3") {
  argx <- x
  freqs <- data.frame(Var1=0:K)
  for (t in type) {
    freq <- data.frame(get(t)(N=N, K=K, m=m))
    names(freq)[2] <- t
    freqs <- merge(freqs, freq, by="Var1", all=T)
  }
  freq <- freqs
  freq[is.na(freq)] <- 0
  freq$binom <- dbinom(0:K, size=K, prob=1/m)
  names(freq)[1] <- "x"
  subset(freq, x %in% argx)
}

#' Theoretical density for lineups under scenario 3 for m = 2, 3, and 20
#' 
#' Some more details to be written later
#' @param x vector of the number of picks of the data plot out of K evaluations
#' @param K number of evaluations
#' @param m size of the lineup
#' @param type one of "mpfr" or "numeric". Should the result be in arbitrary numeric length or be a numeric? Internally the Rmpfr package is used to get maximal precision.
#' @export
#' @examples
#' hdensity(0:5, 5, m=2)
#' ## compare to 
#' dvisual(0:5, 5, m=2)
#' 
#' ## test how many K can be computed without numeric loss
#' for (K in 10:50) { 
#'   print(K); 
#'   print(qplot(0:K, hdensity(0:K, K, m=20))); 
#'   print(sum(hdensity(0:K, K, m=20))); 
#'   scan()
#' }
hdensity <- function(x, K, m, type="numeric") {
  T1 <- function(m) {
    if (m==2) return(expression(1/m*(1/u^2*(log((u+1)/u)*u^2 + u - log(u+1)))))
    if (m==3) return(expression(1/m*(1/u^3*((3*u+1)*log(u) - (u^3-3*u-1)*log((u+1)/u)+ 
                                         (4*u^3-3*u-1)*log((2*u+1)/(2*u)) + 2*u^2 -
                                         3*u*log(2*u) + log(u+1) - log(2*u)))))
  }
  ci <- function(i, K, x) {
    res <- mpfr(rep(0, length=length(i)), 120)
    idx <- which(i >= K-x)
    res[idx] <- ((-1)^(i[idx]-K+x) * chooseMpfr(x, i[idx]-K+x))
#    browser()
    res
  }
  Tone <- function(i, m) {
#    load("data/Dis.RData")
    data(Dis)
    if (i==0) return(1)
    u <- 1
    if (i==1) return(eval(T1(m)))
    (-1)^(i-1)/factorial(i-1)*Dis[[m]][i] #eval(Dks[[i]])
  }
  Tone_saved <- function(i,m) {
    if (i == 0) return(1)
#    load("data/T1m.RData")
    data(T1m, package="vinference")
    T1m[[m]][i]
  }
  Tim <- function(i, m) {
    sapply(i, function(j) Tone_saved(j,m))
  }
  Dk <- function(express, name, k) {
    ## returns all k derivatives of express up to k
    res <- list()
    res[[1]] <- express # zeroeth derivative is expression itself
    for (j in 2:(k+1)) res[[j]] <- D(res[[j-1]], name)
    res
  }
  hone <- function (x, K, m) {
    data(T1m)
    cis <- ci(0:K, K, x)

    #   choose(K, x)*sum(cis*unlist(Tim(0:K, m)))
 #   browser()
    xs <- T1m[[m]][1:K]*chooseMpfr(K, x)*cis[-1]
    sum(xs) + cis[1]
  }
  
  if (!(m %in% c(2,3,20))) {
    stop("Not implemented for this lineup size, use bootstrap simulation with dvisual instead.")  
  }
  if (K > 100) {
    stop("Not implemented for values of K > 100 because of large memory needs. You might want to use bootstrap simulation with dvisual instead.")  
  }
#  Dks <- Dk(T1(m), "u", k=K)
  xs <- lapply(x, hone, K=K, m=m)
  ys <- xs[[1]]
  if (length(xs) > 1)
    for (i in 2:length(xs)) ys <- c(ys, xs[[i]])
  if (type == "mpfr") return(ys)
  else return(as.numeric(ys))
}

#' Quantile of the analytically derived density for visual inference according to scenario 3
#' 
#' @param q (vector) of quantiles
#' @param K number of evaluations
#' @param m lineup size, currently only m=2 and 3 are treated analytically. Use simulation within dvisual to get to other values for m
#' @export
#' @examples
#' ## get critical values of visual triangle test:
#' hquantile(q=c(0.95, 0.99), K=c(5,10,15,20, 25, 30), m=3)
#' 
#' ## get critical values of full lineup test:
#' hquantile(q=c(0.95, 0.99), K=c(5,10,15,20, 25, 30), m=20)
hquantile <- function(q, K, m) {
  dframe <- data.frame(expand.grid(q, K))
  names(dframe) <- c("q", "K")
  require(plyr)
  res <- ddply(dframe, .(q, K), function(x) {
    hs <- cumsum(hdensity(x=0:x$K, K=x$K, m=m))
    which(hs>=x$q)[1]
  })
  names(res)[3] <- "x"
  res$x[res$x > res$K] <- NA
  res
}
  
#' Explicit density function of visual inference under scenario 2 for m = 3
#'
#' more details to follow
#' @param x number of times data plot was picked
#' @param K number of  evaluations by independent observers
#' @param m lineup size. Only implemented for m=3
#' @export
dv2 <- function(x,K, m=3) { 
  dv2one <- function(x, K, m=m) {
    g <- function(q) {
      res <- q*((q+2)*log(q+2)-2*(q+1)*log(q+1)+q*log(q))
      res[q==0] <- 0
      res
    }
    if (m!=3) stop("Function not implemented for lineup sizes <> 3. Use bootstrap based density estimation in dvisual instead.")
    
    f <- function(q, K, x) choose(K,x)*g(q)^x*(1-g(q))^(K-x)
    integrate(f, 0,1,K=K,x=x)$value
  }
  unlist(sapply(x, dv2one, K=K, m=m))
}

#' Explicit density function of visual inference under scenario 3 for m = 2
#' 
#' @param x
#' @param K
#' @export
#' @examples
#' h(0:5, K = 5)
#' ## compare to 
#' dvisual(0:5, K=5, m=2)
h <- function(x, K) {
  hone <- function(x, K) {
    S <- function(m, n, a, b) {
 #     print(match.call())
      C <- function(m, n, i) {
  #      print(match.call())
        if (i >= m) return(1)
        j <- (i+1):m
        
   #     print(prod(j/(j+n-m-1)))
        return(prod(j/(j+n-m-1)))
      }
      D <- function(m,n,a,b) {
    #    print(match.call())
        db <- b^m/(1+b)^{n-1}
        da <- a^m/(1+a)^{n-1}
        if (b == Inf) {
          db <- 0
        }
     #   print(1/(n-1)*(db - da))  
        return (1/(n-1)*(db - da))
      }
      if (m == 0) {
        if (n == 0) return (b-a)
        if (n == 1) return (log((1+b)/(1+a)))
        return(-D(0, n, a, b))
      }
      if (m == 1) {
        if (n == 1) return(b-a - log((1+b)/(1+a)))
      }
      
      if (m == n) {
        sumD <- 0
        for (i in 2:m) sumD <- sumD + D(i,i, a, b)*C(m, m, i)
        return(C(m,m,1)*S(1,1,a,b) - sumD)
      }
    
      sumD <- 0
      for (i in 1:m) sumD <- sumD + D(i,n+i-m, a, b)*C(m, n, i)
      return(C(m,n,0)*S(0,n-m,a,b) - sumD)      
    #  return(-D(m,n,a,b)+m/(n-1)*S(m-1, n-1, a, b))
    }
    if (x < 2) x <- K-x
    choose(K, x) *(0.5*S(x, K, 0, 1) + 0.5*S(x-2, K, 1, Inf))
  }
  sapply(x, hone, K=K)
}


#' List of coefficients to evaluate hdensity
#' 
#' @name T1m
#' @title List of coefficients in hdensity
#' @description List of theoretical coefficients to evaluate hdensity in cases m= 2, 3 and 20
#' @docType data
#' @usage data(T1m)
NULL

#' Theoretical distribution for lineups under scenario 3 for m = 2, 3, and 20
#' 
#' Some more details to be written later
#' @param x vector of the number of picks of the data plot out of K evaluations
#' @param K number of evaluations
#' @param m size of the lineup
#' @param type one of "mpfr" or "numeric". Should the result be in arbitrary numeric length or be a numeric? Internally the Rmpfr package is used to get maximal precision.
#' @export
#' @examples
#' hdistribution(0:5, 5, m=3)
hdistribution <- function(x, K, m, type="numeric") {
  hdone <- function(x1, K, m) {
    sum(hdensity(x1:K, K, m))
  }
  sapply(x, hdone, K=K, m=m)
}