#' Calculate pairwise co-membership probability
#'
#' Given a vector of cluster numbers and projection clustering output, \code{postPairs}
#' calculate the posterior probability of any pair of subjects being clustered in the
#' same group.
#'
#' @inheritParams clustKL
#' @inheritParams postMean
#' @param nDraw Number of draws to sample from projection clustering output
#' @param nClusters Vector of cluster numbers
#'
#' @return List of 2 items.
#' ls_prob: pairwise probability tables corresponding to random effect projection
#' specified in \code{ls_idxA}. Row and column index of table indicate subject
#' ID number in \code{dat}. Only the lower triangular matrix is filled.
#' arr_cluster: array of optimized cluster labels based on randomly drawn samples
#' corresponding to chosen random effects.
#' @export
#'
#' @examples
#' data(df_of_draws)
#' ls_idxA <- list(
#'   seq(10),
#'   1:4,
#'   5:7,
#'   8:10
#' )
#' out_pc <- postPairs(df_of_draws,
#'   x_var = paste0("Z", 1:10), id_var = "ID", dat = DATASET,
#'   ls_idxA, nIter = 10, nDraw = 2, nClusters = 4, regQ = 1e-6, seed = 1
#' )
postPairs <- function(df_of_draws, x_var, id_var = "ID", dat, ls_idxA, nIter = 10, nDraw = 1000, nClusters, regQ = 1e-6, seed = 1) {
  id <- as.numeric(dat[, id_var])
  nSub <- length(unique(id))
  nBasis <- length(x_var)

  if (length(nClusters) == 1) nClusters <- rep(nClusters, length(ls_idxA))

  foreach(t = seq(1, length(ls_idxA))) %dopar% {
    idxA <- ls_idxA[[t]]
    idxB <- setdiff(ls_idxA[[1]], idxA)
    nCluster <- nClusters[t]

    set.seed(seed)

    cat(sprintf("Set %d: drawing samples from posterior\n", t))
    post.obj <-
      post.cluster(df_of_draws, nDraw, dat, nBasis, idxA, idxB, id, nCluster, nIter, regQ, seed)

    # posterior probability of 2 subjects being in 1 cluster
    tb <- matrix(NA, nrow = nSub, ncol = nSub)
    for (i in seq(2, nSub)) {
      for (j in seq(i - 1)) {
        tb[i, j] <- mean(post.obj$clusterMatrix[i, ] == post.obj$clusterMatrix[j, ])
      }
    }
    ls_prob_t <- tb
    list(
      ls_prob_t,
      post.obj$clusterVec
    )
  } -> temp

  list(
    ls_prob = lapply(temp, "[[", 1),
    arr_cluster = array(unlist(lapply(temp, "[[", 2)),
      dim = c(nSub, length(ls_idxA))
    )
  )
}

# main helper function
post.cluster <- function(df_of_draws, nDraw, dat, nBasis, idxA, idxB, id, nCluster, nIter, regQ, seed) {
  nSub <- length(unique(id))
  idxDraw <- sample(nrow(df_of_draws), nDraw)
  idx <- 0
  foreach(idx = idxDraw) %dopar% {
    bAMatrix <- matrix(unlist(
      df_of_draws[idx, sprintf("beta[%d,%d]", rep(seq(nSub), each = nBasis), seq(nBasis))]
    ),
    nrow = nSub, byrow = T
    )[, idxA]
    Gamma <- unlist(df_of_draws[idx, "sigmaepsilon"])^2
    G <- matrix(unlist(
      df_of_draws[idx, sprintf("VarBeta[%d,%d]", rep(seq(nBasis), each = nBasis), seq(nBasis))]
    ),
    nrow = nBasis, byrow = T
    )
    Q <- Q.func(idxA, idxB, dat, Gamma, G, id)

    temp <- cluster.iter(NULL, nIter, nCluster, Q, bAMatrix, regQ, seed, do.plot = F)
    dK <- temp[[2]][, , nIter]
    clusterVec <- temp[[1]][, nIter]

    list(KLK.func(dK, clusterVec, bAMatrix, Q), clusterVec, dK, Q, bAMatrix)
  } -> temp

  ls_Q <- array(unlist(lapply(temp, "[[", 4)), dim = c(nSub, length(idxA), length(idxA), nDraw))
  ls_bAMatrix <- array(unlist(lapply(temp, "[[", 5)), dim = c(nSub, length(idxA), nDraw))
  temp1 <- cluster.iter(NULL, nIter, nCluster, ls_Q, ls_bAMatrix, regQ, seed, do.plot = F)
  clusterVec <- temp1[[1]][, nIter]

  list( # must check if element sequence is right
    KL = unlist(lapply(temp, "[[", 1)),
    clusterMatrix = matrix(unlist(lapply(temp, "[[", 2)), ncol = nDraw, byrow = F),
    dKArray = array(unlist(lapply(temp, "[[", 3)), dim = c(nCluster, length(idxA), nDraw)),
    clusterVec = clusterVec
  )
}
