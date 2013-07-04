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
  table(replicate(N,  {
    dataprob <- runif(1)
    nulls <- runif(K-1)
    
    individual <- sum(replicate(K, lineup(m, dataprob=dataprob, nulls=nulls))==1)
    individual
  }))/N
}


scenario4 <- function(N, K, m=20) {
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
#' @param type type of simulation used: scenario 4 assumes that the same lineup is shown in all K evaluations
#' @param upper.tail compute probabilities P(X >= x). Be aware that the use of this parameter is not consistent with the other distribution functions in base. There, a value of P(X > x) is computed for upper.tail=TRUE.
#' @return simulation based p value to observe x or more picks of the data plot in K evaluation under the assumption that the data plot is consistent with the null hypothesis. For comparison a p value based on a binomial distribution is provided as well.
#' @export
#' @examples
#' pvisual(15, 20, m=3) # triangle test
pvisual <- function(x, K, m=20, N=10000, type="scenario4", upper.tail=TRUE) {
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
#' @param type type of simulation used: scenario 4 assumes that the same lineup is shown in all K evaluations
#' @return simulation based density to observe x picks of the data plot in K evaluation under the assumption that the data plot is consistent with the null hypothesis. For comparison a p value based on a binomial distribution is provided as well.
#' @export
#' @examples
#' dvisual(2, 20, m=3) # triangle test
dvisual <- function(x, K, m=20, N=10000, type="scenario4") {
  freq <- get(type)(N=N, K=K, m=m)
  sim <- sapply(x, function(y) freq[as.numeric(names(freq)) %in% y])
  return(cbind(x=x, "simulated"=sim, "binom"=dbinom(x, size=K, prob=1/m)))
}
