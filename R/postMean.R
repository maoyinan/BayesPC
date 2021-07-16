#' Calculate posterior mean from the Linear Mixed Model result
#'
#' Compute parameters' posterior mean based on \link{modelStan} output.
#' @param df_of_draws Data frame of simulated LMM output
#' @inheritParams modelStan
#' @return list(b0,bMatrix,G,Gamma)
#'
#' @family BayesPC main functions
#' @export
#'
#' @examples
#' \dontrun{
#' df_of_draws <- modelStan("Record", paste0("Z", 1:10), "ID", DATASET)
#' }
#' ls_par <- postMean(df_of_draws, paste0("Z", 1:10), "ID", DATASET)
postMean <- function(df_of_draws, x_var, id_var, dat) {
  nBasis <- length(x_var)
  nSub <- length(unique(dat[, id_var]))
  b0 <- mean(df_of_draws[, "beta0"])
  bMatrix <- matrix(
    apply(
      df_of_draws[, sprintf("beta[%d,%d]", rep(seq(nSub), each = nBasis), seq(nBasis))],
      2, mean
    ),
    nrow = nSub, byrow = T
  )
  G <- matrix(
    apply(
      df_of_draws[, sprintf("VarBeta[%d,%d]", rep(seq(nBasis), each = nBasis), seq(nBasis))],
      2, mean
    ),
    nrow = nBasis, byrow = T
  )
  Gamma <-
    mean(df_of_draws[, "sigmaepsilon"])^2

  list(b0 = b0, bMatrix = bMatrix, G = G, Gamma = Gamma)
}
