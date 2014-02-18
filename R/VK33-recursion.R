p3 <- function(K, x) {
  require(plyr)
  as.vector(choose(K,x)*unlist(llply(x, function(y) h3(x=y, K=K))))
}

h3 <- function(K, x) {
  if (K < 0) return(NA)
  if (K==1) {
    if (x == 0) return(2/3)
    if (x == 1) return(1/3)
    else return(NA)
  }
  if (x == 0) {
    return(( 1- g3(K, 0))/(K-1) )
  }
  return((x*h3(K-1, x-1) - g3(K, x))/(K-1))
}



  
g3 <- function(K, x) {
  if (K < 2) return(NA)
  if (K == 2) {
    if (x == 0) return(log(27/16))
    if (x == 1) return(1 - log(27/16))
    if (x == 2) return(log(27/16))
    else return(NA)
  }
  if (x == K) {
    if(K==3) return(f3(3,3))
    return((27/3^K - 16/2^K + 1)/((K-2)*(K-3)))
  }
  return((f3(K, x) + (K-x)*g3(K-1, x))/(K-2))
}

f3 <- function(K, x) {
  if (K < 2) return (NA)
  if (K == 2) {
    if (x == 0) return(-2)
    if (x == 1) return(-1)
    if (x == 2) return (0)
    else return(NA)
  }
  if (K==3) {
    if (x==0) return(-log(4/3)-1)
    if (x==1) return (log(4/3) - 1)
    if (x==2) return(-log(4/3))
    if (x==3) return(log(4/3))
    else return(NA)
  }
  if (x == K) {
    return( (1-2/2^(K-3) + 1/3^(K-3))/(K-3))
  }
  return((-2^(4-K) + 2^(K-x)/3^(K-3) + (K-x)*f3(K-1, x))/(K-3))
}