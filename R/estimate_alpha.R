#' Functions to produce visual estimate of alpha


#' Simulate data from a lineup evaluation experiment using the Dirichlet-Multinomial model
#' 
#' @param alpha The Dirichlet parameter which is related to the number of interesting panels
#' @param m0 The number of null panels in the lineup
#' @param K The total number of null panel selections (or, in a Rorschach lineup, the total number of evaluations)
#' @param N Number of lineups to simulate
#' @importFrom gtools rdirichlet
#' @importFrom stats rmultinom
sim_lineup_model <- function(alpha, m0 = 19, K = 22, N = 50) {
  theta <- gtools::rdirichlet(1, rep(alpha, m0))
  sels <- stats::rmultinom(N, size = K, prob = theta)
  sels
}

#' Compute the expected number of c-interesting panels for a lineup experiment
#' 
#' @param alpha The Dirichlet parameter which is related to the number of interesting panels
#' @param c The number of selections a panel must have to be interesting (can be non-integer)
#' @param m0 The number of null panels in the lineup
#' @param K The total number of null panel selections (or, in a Rorschach lineup, the total number of evaluations)
#' @export
expected_number_panels <- function(alpha, c=m0/K, m0 = 19, K=30) {
  x <- ceiling(c):K
  summation <- choose(K, x) * beta(x + alpha, K - x + (m0 - 1)*alpha)
  
  m0/beta(alpha, (m0 - 1)*alpha)*sum(summation)
}

#' Simulate the number of c-interesting panels for a lineup experiment
#' 
#' @param alphas Numeric vector of alpha values to conduct simulations for
#' @param c The number of selections a panel must have to be interesting (can be non-integer)
#' @param m0 The number of null panels in the lineup
#' @param K The total number of null panel selections (or, in a Rorschach lineup, the total number of evaluations)
#' @param N_points The number of points to simulate for each value of alpha
#' @param avg_n_sims The number of simulations to average to get a single point value.
#'          Averaging several simulations reduces the visual noise but also decreases 
#'          the separation between possible values for a more continuous appearance.
#' @export
#' @importFrom tidyr unnest
#' @importFrom tibble tibble
#' @importFrom purrr map
#' @importFrom dplyr mutate group_by summarize
#' @examples 
#' sim_interesting_panels()
sim_interesting_panels <- function(alphas = 10^seq(-2, 2, .05), c = m0/K, m0 = 19, K = 30, 
                                   N_points = 10, avg_n_sims = 10) {
  # Each point (of N_points) is an average of avg_n_sims separate lineups
  # First, generate all of the lineups and count the number of interesting panels
  df <- tibble::tibble(
    alpha = alphas,
    plot_sels = purrr::map(.data$alpha, sim_lineup_model, N = N_points*avg_n_sims, K = K),
    interesting_panels = purrr::map(.data$plot_sels,
      ~tibble::tibble(n_interesting = colSums(.x >= c),
                      rep = 1:ncol(.x) - 1))
  ) %>% 
    tidyr::unnest(.data$interesting_panels)
  
  # Then, average the panels together a bit
  df2 <- df %>%
    dplyr::mutate(point_num = (rep - (rep %% avg_n_sims))/N_points) %>%
    dplyr::group_by(.data$alpha, .data$point_num) %>%
    dplyr::summarize(n_interesting = mean(.data$n_interesting))
  
  df2
}

#' Create a plot for visual estimation of alpha
#' 
#' The visual estimation plot contains the expected number of interesting 
#' panels for each alpha value under the specified experimental conditions and
#' simulated lineup experiments under each scenario which provide some 
#' information about the expected variability. 
#' 
#' @param c The number of selections a panel must have to be interesting (can be non-integer)
#' @param m0 The number of null panels in the lineup
#' @param K The total number of null panel selections (or, in a Rorschach lineup, the total number of evaluations)
#' @param alphas Numeric vector of alpha values to conduct simulations for
#' @param ... additional arguments to sim_interesting_panels
#' @export
#' @importFrom dplyr arrange mutate `%>%` count group_by
#' @importFrom purrr map_dbl
#' @importFrom tibble tibble
#' @importFrom tidyr crossing
#' @import ggplot2
#' @examples 
#' estimate_alpha_visual()
estimate_alpha_visual <- function(c = m0/K, m0 = 19, K = 30, alphas = 10^seq(-3, 2, .05), ...) {
  n <- z <- NULL
  
  # Get theoretical function
  model_df <- tibble::tibble(alpha = alphas,
                             n_sel_plots = alphas %>% 
                               purrr::map_dbl(expected_number_panels, c = c, m0=m0, K=K)
  )
  
  # Get simulated data
  prior_pred_mean <- sim_interesting_panels(c = c, m0= m0, K = K, alphas = alphas, ...) %>%
    dplyr::arrange(.data$alpha) %>%
    dplyr::mutate(label = sprintf("alpha == %f", .data$alpha) %>%
                    factor(levels = sprintf("alpha == %f", alphas), ordered = T))
  
  # Compute # points from prior_pred_mean
  pts_per_alpha <- prior_pred_mean %>% dplyr::group_by(.data$alpha) %>% dplyr::count() %>% `[[`("n") %>% mean
  
  # minor break points
  mb <- tidyr::crossing(x = c(2.5, 5, 7.5), y = 10^(seq(-4, 5))) %>%
    dplyr::mutate(z = .data$x*.data$y) %>%
    `[`("z") %>%
    unlist() %>%
    as.numeric
  
  bks <- 10^seq(-4, 5, 1)
  
  # Finally, the actual plot
  ggplot() +
    geom_point(aes(x = .data$alpha, y = .data$n_interesting),data = prior_pred_mean, alpha = 1/pts_per_alpha) +
    geom_line(aes(x = .data$alpha, y = .data$n_sel_plots), data = model_df, color = "blue", size = 1) +
    scale_x_log10(name = expression(alpha), breaks = bks,
                  minor_breaks = mb,
                  labels = bks) +
    coord_cartesian(xlim = range(alphas)) +
    scale_y_continuous(name = sprintf("Average number of panels with at least %.2f selection(s)", c), breaks = 1:K)
  
}


#' Numerically estimate alpha using the average number of c-interesting panels
#' 
#' This function calculates a precise numerical estimate of alpha based by matching 
#' the number of observed interesting panels to the corresponding value of alpha.
#' 
#' @param Zc Average number of panels with at least c selections
#' @param c The number of selections a panel must have to be interesting (can be non-integer)
#' @param m0 The number of null panels in the lineup
#' @param K The total number of null panel selections (or, in a Rorschach lineup, the total number of evaluations)
#' @return a data frame with list-columns for alpha, optimized sum of squared error, counts (as returned from optim), convergence, 
#' @export
#' @importFrom purrr map map_dbl map_int map_chr
#' @importFrom tibble tibble
#' @examples 
#' res <- estimate_alpha_numeric(5, c = 1, m0 = 19, K = 30)
#' res$alpha
estimate_alpha_numeric <- function(Zc, c = m0/K, m0 = 19, K = 30) {
  # Just a wrapper
  purrr::map(Zc, estimate_alpha_num_scalar, c = c, m0 = m0, K = K) %>% {
    tibble(
      alpha = purrr::map_dbl(., "alpha"),
      sum_sq_error = purrr::map_dbl(., "sum_sq_error"),
      counts = purrr::map(., "counts"),
      convergence = purrr::map_int(., "convergence"),
      message = purrr::map_chr(., ~ifelse(length(.$message) == 0, NA, .$message))
    )}
}

#' @importFrom tibble tibble
#' @importFrom dplyr mutate filter
#' @importFrom stats optim
estimate_alpha_num_scalar <- function(Zc, c = m0/K, m0 = 19, K = 30) {
  stopifnot(Zc < K, Zc > 0)
  stopifnot(c > 0, m0 > 1, K > 1)
  
  
  inv_exp_panels <- function(alpha, X, c, m0, K) {
    if (alpha <= 0) return(Inf)
    
    (X - expected_number_panels(alpha = alpha, c = c, m0 = m0, K = K))^2
  }
  
  # Get good initialization values
  df <- tibble(alpha = 10^seq(-2, 2, .5)) %>%
    dplyr::mutate(objval = purrr::map_dbl(.data$alpha, inv_exp_panels, X = Zc, c = c, m0 = m0, K = K)) %>%
    dplyr::filter(.data$objval == min(.data$objval))

  
  res <- optim(list(alpha = df$alpha), inv_exp_panels, 
               X = Zc, c = c, m0 = m0, K = K, 
               method = "Brent", lower = 1e-4, upper = 100)
  
  names(res) <- c("alpha", "sum_sq_error", "counts", "convergence", "message")
  
  if (res$alpha < 0.01) warning("Warning: alpha estimate is too low to be reliable. Null panel generation method may produce null plots which are too visually distinct.")
  
  res
}

#' Add a band corresponding to the observed number of interesting panels to the visual estimation plot
#' 
#' This function adds a band to assist in the visual estimation of alpha from evaluated lineups. 
#' 
#' @param obs Observed number of "c-interesting" panels (panels with at least c selections). If NA, the line is omitted.
#' @param limits Lower and upper bounds for number of panels. If NA, the bands are omitted.
#' @param c The number of selections a panel must have to be interesting (can be non-integer)
#' @param m0 The number of null panels in the lineup
#' @param K The total number of null panel selections (or, in a Rorschach lineup, the total number of evaluations)
#' @export
#' @importFrom tibble tibble
#' @importFrom ggplot2 geom_line geom_polygon
#' @importFrom purrr map_dbl
#' @examples
#' estimate_alpha_visual() + observed_band(obs = NA, limits = c(5, 7))
#' estimate_alpha_visual() + observed_band(obs = 6, limits = NA)
#' estimate_alpha_visual(c = 1, m0 = 19, K = 50) + observed_band(obs = 6, limits = c(4, 8), c = 1, m0 = 19, K = 50)
observed_band <- function(obs, limits, c = m0/K, m0 = 19, K = 30) {

  # get alpha values for critical y values
  obs_line <- if (!is.na(obs)) {
    obs_alph <- estimate_alpha_num_scalar(obs, c = c, m0 = m0, K = K)$alpha
    
    # Construct lines for center/observed
    segs <- rbind(
      tibble(x = c(1/Inf, obs_alph[1]), y = obs, type = "horiz"),
      tibble(x = obs_alph[1], y = c(obs, -Inf), type = "vert")
    )
    
    geom_line(aes(x = .data$x, y = .data$y, group = .data$type), data = segs, color = "grey30", size = 0.5)
  } else {
    NULL
  }
  
  obs_band <- if (!any(is.na(limits))) {
    limits_alph <- estimate_alpha_numeric(limits, c = c, m0 = m0, K = K)$alpha
    
    # Construct polygon coords for shading
    bands <- rbind(
      tibble(x = c(1/Inf, limits_alph, 1/Inf, 1/Inf),
             y = limits[c(1, 1, 2, 2, 1)],
             type = "horiz"),
      tibble(x = limits_alph[c(1, 1, 2, 2, 1)],
             y = c(-Inf, limits, -Inf, -Inf),
             type = "vert")
    )
    
    geom_polygon(aes(x = .data$x, y = .data$y, group = .data$type), data = bands, fill = "grey", alpha = 0.5)
  } else {
    NULL
  }
  return(list(obs_line, obs_band))
}

