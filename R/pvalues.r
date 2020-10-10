#' Visual Inference for different Lineup Scenarios under a Dirichlet Model
#' 
#' Density, distribution function and quantiles for visual inference scenarios.
#' Visual inference is used to determine significance of a visual finding. 
#' The lineup protocol (Buja et al., 2009) establishes a formal framework for 
#' testing graphical findings. The package \code{nullabor} helps with the creation 
#' of lineups using various null generation models.
#' Here, we provide functions to evaluate results.
#' 
#' When administering visual tests, we distinguish between three different scenarios:
#' \itemize{
#' \item Scenario 1: in each of K evaluations a different data set and a different set of (m-1) null plots is shown.
#' \item Scenario 2: in each of K evaluations the same data set but a different set of (m-1) null plots is shown.
#' \item Scenario 3: the same lineup, i.e. same data and same set of null plots, is shown to K different observers.
#' }
#' Under scenario 3, the number of data picks under the null hypothesis that the data plot is 
#' visually not more salient than one of the null plots is distributed according 
#' to a ratio of Beta functions:
#' \deqn{P (X = x) = {K \choose x} \frac{B(x + \alpha, K-x+(m-1)\alpha)}{B(\alpha, (m-1)\alpha)}}
#' where 
#' \eqn{B(.,.)} is the Beta function, 
#' \eqn{\alpha > 0} is the rate of a flat Dirichlet distribution,
#' \eqn{K} is the number of times the lineup has been evaluated, 
#' \eqn{x} number of times the data plot has been picked as the visually most interesting,
#' \eqn{m} is the number of panels in a lineup (the lineup size).
#' 
#' For large values of alpha, scenario 3 converges to scenario 1. 
#' @param x vector, number of data identifications,
#' @param p (vector of) probabilities,
#' @param K positive value, number of evaluations of the lineup,
#' @param m number of panels in the lineup,
#' @param alpha positive value, rate parameter of the flat Dirichlet distribution,
#' @param scenario integer value.
#' @param lower.tail defaults to TRUE, if TRUE probabilities are \eqn{P(X ≤ x)}, otherwise, \eqn{P(X ≥ x)}. 
#' Note that the second probability is a deviation from R standard: usually \eqn{P(X > x)} is returned. 
#' However, here, returning \eqn{P(X ≥ x)} is more useful in an inference setting, as it corresponds to the 
#' p value.
#' @importFrom stats dbinom pbinom qbinom
#' @return The functions return (vector) of probabilities for \eqn{P(X = x)}, \eqn{P(X ≤ x)}, or quantiles.  
#' @export
#' @references 
#' Andreas Buja, Dianne Cook, Heike Hofmann, Michael Lawrence, Eun-Kyung Lee, Deborah F. Swayne and Hadley Wickham, Statistical inference for exploratory data analysis and model diagnostics.
#' Phil. Trans. R. Soc. A. 367: 4361-4383, 2009, \url{https://doi.org/10.1098/rsta.2009.0120}
#' @examples 
#' # Probabilities to see between 5 and 10 data identifications
#' # in lineup of size 20, with 15 evaluations, and an estimated 
#' # alpha of 0.1
#' dVis(x=5:10, K= 15, m=20, alpha = 0.1)
#' 
#' dframe <- data.frame(
#'   x=0:15, 
#'   probabilities = c(dVis(0:15, K=15, alpha = 0.1),
#'                     dVis(0:15, K=15, alpha = 1),
#'                     dVis(0:15, K=15, alpha = 5)),
#'   alpha = factor(rep(c(0.1, 1, 5), each = 16))                  
#' )
#' library(ggplot2)
#' ggplot(data = dframe, aes(x = x, y = probabilities, colour = alpha)) +
#'   geom_point() 
#'   
#'  # how many data picks do we need in a lineup of size m
#'  # with K = 30 evaluations and alpha = 0.1 to achieve a 
#'  # significance at 10%, 5% or 1%?
#'  qVis(p = c(0.9, 0.95, 0.99), K = 30, m=20, alpha = 0.1)
dVis <- function(x, K, m = 20, alpha, scenario = 3) {
  if (scenario == 3) {
    return(choose(K, x) * beta(x+alpha, K-x+(m-1)*alpha)/beta(alpha, (m-1)*alpha))
  }
  if (scenario == 1) { # binomial distribution
    return(dbinom(x, K, prob = 1/m))
  }
}


#' @rdname dVis
#' @export
pVis <- function(x, K, m = 20, alpha, scenario = 3, lower.tail = TRUE) {
  if (scenario == 1) { # binomial distribution
    if (lower.tail == FALSE) x <- x-1
    return(pbinom(x, size = K, prob = 1/m, lower.tail = lower.tail))
  }
  
  if (scenario == 3) {
    # x should be single value now, or we need to be a bit inventive
    if (lower.tail) {
      xs <- seq(0, max(x, na.rm=TRUE))
    } else (
      xs <- seq(K:min(x, na.rm=TRUE))
    )
    dvis <- dVis(x = xs, K =K, m=m, alpha = alpha, scenario=scenario)
    cum <- cumsum(dvis)
    if (!lower.tail) cum <- rev(cum)
      
    return(cum[x+1])
  }
}

#' @rdname dVis
#' @export
qVis <- function(p, K, m = 20, alpha, scenario = 3) {
  if (scenario == 1) { # binomial distribution
    return(qbinom(p=p, size = K, prob = 1/m))
  }
  
  if (scenario == 3) {
    xs <- seq(0, K)
    # x should be single value now, or we need to be a bit inventive
    
    dvis <- dVis(x = xs, K =K, m=m, alpha = alpha, scenario=scenario)
    cum <- cumsum(dvis)
    idx <- sapply(p, FUN=function(pp) {
        wh <- which(pp >= cum)
        if (length(wh)!=0) max(wh)
        else 0
      })
    
    return(idx)
  }
}



