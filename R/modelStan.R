#' Bayesian mixed linear regression modelling
#'
#' To fit longitudinal dataset using linear mixed model by calling stan package.
#'
#' @param y_var character of response variable
#' @param x_var character vector of random effect variables
#' @param id_var character of id variable
#' @param dat dataset to model with stan
#' @param ... Additional arguments passed to \link[rstan]{sampling}
#'
#' @return dataframe of posterior output from Bayesian mixed linear regression
#' @family BayesPC main functions
#' @import rstan
#' @export
#' @examples
#' data(DATASET)
#' df_of_draws <- modelStan("Record", paste0("Z", 1:10), "ID", DATASET, iter = 2, chain = 1)
modelStan <- function(y_var, x_var, id_var, dat, ...) {
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
  system.time(fit02 <- sampling(sm, data = dat_mcmc, seed = 1, cores = 1, init = 0, ...))
  df_of_draws <- as.data.frame(fit02)
  df_of_draws
}
