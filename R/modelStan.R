#' Bayesian mixed linear regression modelling
#'
#' To fit longitudinal dataset using linear mixed model by calling stan package.
#'
#' @param y_var Character of response variable
#' @param x_var Character vector of random effect variables
#' @param id_var Character of id variable
#' @param dat Longitudinal data input
#' @param seed Random seed to pass to \link[rstan]{sampling}
#' @param ... Additional arguments passed to \link[rstan]{sampling}
#'
#' @return dataframe of posterior output from Bayesian mixed linear regression
#' @family BayesPC main functions
#' @import rstan
#' @export
#' @examples
#' df_of_draws <- modelStan("Record", paste0("Z", 1:10), "ID", DATASET, seed=1, iter = 2, chain = 1)
modelStan <- function(y_var, x_var, id_var="ID", dat, seed=1, ...) {
  dat_mcmc <- list(
    N = nrow(dat),
    y = dat[, y_var],
    nBasis = length(x_var),
    Z = dat[, x_var],
    nSub = length(unique(dat[, id_var])),
    sub = as.numeric(dat[, id_var])
  )
  rt <- stanc(file = system.file("extdata", "mcmc_try2.stan", package = "BayesPC"))
  sm <- stan_model(stanc_ret = rt, verbose = FALSE)
  system.time(fit02 <- sampling(sm, data = dat_mcmc, seed = seed, cores = 1, init = 0, ...))
  df_of_draws <- as.data.frame(fit02)
  df_of_draws
}
