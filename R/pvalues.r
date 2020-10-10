#' Visual Inference
#' 
#' XXX intro sentence for Visual inference
#' We distinguish between three different scenarios:
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
#' For large values of alpha, scenario 3 converges to scenario 1. 
#' @param x vector, number of data identifications,
#' @param K positive value, number of evaluations of the lineup,
#' @param m number of panels in the lineup,
#' @param alpha positive value, rate parameter of the flat Dirichlet distribution,
#' @param scenario integer value.
#' @importFrom stats dbinom
#' @export
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
#' dframe %>% 
#'   ggplot(aes(x = x, y = probabilities, colour = alpha)) +
#'   geom_point() 
dVis <- function(x, K, m = 20, alpha, scenario = 3) {
  if (scenario == 3) {
    choose(K, x) * beta(x+alpha, K-x+(m-1)*alpha)/beta(alpha, (m-1)*alpha)
  }
  if (scenario == 1) { # binomial distribution
    dbinom(x, K, prob = 1/m)
  }
}