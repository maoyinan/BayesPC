#' Projection clustering function
#'
#' Compute LMM fitted values based on posterior mean output of parameters,
#' with chosen random effects to project on.
#'
#' @param ls_par Posterior mean of LMM parameters generated from \link{postMean} output
#' @inheritParams modelStan
#' @param ls_idxA List of random effect indices to project on
#' @family BayesPC main functions
#'
#' @return Matrix containing fitted values, with each column
#' corresponding to a set of chosen random effects indices
#' specified in \code{ls_idxA}
#' @import foreach
#' @export
#'
#' @examples
#' \dontrun{
#' df_of_draws <- modelStan("Record", paste0("Z", 1:10), "ID", DATASET)
#' }
#' ls_par <- postMean(df_of_draws, paste0("Z", 1:10), "ID", DATASET)
#' ls_idxA <- list(
#'   seq(10),
#'   1:4,
#'   5:7,
#'   8:10
#' )
#' mat_fitted <- pcFit("ID", ls_par, DATASET, ls_idxA)
pcFit <- function(id_var, ls_par, dat, ls_idxA) {
  mat_fitted <- data.frame(
    matrix(NA, nrow = nrow(dat), ncol = length(ls_idxA), dimnames = list(c(), sprintf("y_hatA%d", seq(length(ls_idxA)))))
  )
  id <- as.numeric(dat[, id_var])

  b0 <- ls_par$b0
  bMatrix <- ls_par$bMatrix
  G <- ls_par$G
  Gamma <- ls_par$Gamma


  foreach(t = seq(1, length(ls_idxA))) %dopar% {
    idxA <- ls_idxA[[t]]
    idxB <- setdiff(ls_idxA[[1]], idxA)

    bAMatrix <- bMatrix[, idxA]
    y <- my.fit(idxA, idxB, dat, G, bAMatrix, b0, id)
    list(y)
  } -> temp
  mat_fitted <- matrix(unlist(lapply(temp, "[[", 1)), ncol = length(ls_idxA), byrow = F, dimnames = list(c(), sprintf("y_hatA%d", seq(length(ls_idxA)))))
}

# helper function
my.fit <- function(idxA, idxB, dat, G, bAMatrix, b0, id) {
  nSub <- length(unique(id))
  GA <- G[idxA, idxA]
  GA1 <- solve(GA)
  if (length(idxB) > 0) {
    GAB <- G[idxA, idxB]
    ZB <- as.matrix(dat[, paste0("Z", idxB)])
    term <- t(GAB) %*% GA1
  }
  ZA <- as.matrix(dat[, paste0("Z", idxA)])

  y <- numeric(length(id))
  for (i in seq(length(y))) {
    if (length(idxB) > 0) {
      y[i] <- b0 + as.numeric((ZA[i, ] + ZB[i, ] %*% term) %*% bAMatrix[id[i], ])
    } else {
      y[i] <- b0 + as.numeric((ZA[i, ]) %*% bAMatrix[id[i], ])
    }
  }
  y
}
