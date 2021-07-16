#' Choose cluster number by KL method
#'
#' To find minimal cluster number which reduces Kullback-Leibler divergence
#' no larger than a customized threshold. KL monotonically decreases with number
#' of clusters.
#'
#' @inheritParams pcFit
#' @param nIter Number of iterations used in clustering optimization
#' @param thKL Proportion of max KL, used as upper threshold of Kullback-Leibler divergence corresponding to chosen number of clusters
#' @param regQ Positive regularization value to add to the diagonal of matrix to be inverted
#' @param seed Random seed to initialize cluster centers
#'
#' @return list(nClusters,KLs,cluster0)
#' \item{nClusters}{Vector of optimised cluster numbers corresponding to sets of random effects in \code{ls_idxA}}
#' \item{KLs}{List of Kullback-Leibler divergence across grid of cluster numbers corresponding to \code{ls_idxA}}
#' \item{cluster0}{List of matrix, column k indicates cluster membership for k clusters corresponding to \code{ls_idxA}}
#'
#' @family BayesPC main functions
#' @family Cluster number choices
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
#' out_KL <- clustKL("ID", ls_par, DATASET, ls_idxA, 10, .1, 1e-6, 1)
clustKL <- function(id_var, ls_par, dat, ls_idxA, nIter, thKL, regQ, seed){
  id <- as.numeric(dat[, id_var])
  nSub <- length(unique(id))

  bMatrix <- ls_par$bMatrix
  G <- ls_par$G
  Gamma <- ls_par$Gamma

  foreach(t = seq(1,length(ls_idxA))) %dopar% {
    idxA <- ls_idxA[[t]]
    idxB <- setdiff(ls_idxA[[1]],idxA)

    cat(sprintf('Set %d: optimizing cluster number by KL\n',t))

    bAMatrix <- bMatrix[,idxA]
    Q <- Q.func(idxA, idxB, dat,Gamma, G,id)
    Kcluster.obj <- KL.pick(nIter, bAMatrix,Q, regQ, do.plot=F, seed=seed)
    threshold <- Kcluster.obj$KL[1]*thKL
    nClusters_t <- which(Kcluster.obj$KL<threshold)[1]
    KLs_t <- Kcluster.obj$KL
    cluster_t <- Kcluster.obj$clusterMatrix

    list(nClusters_t, KLs_t, cluster_t)
  }->temp

  list(
    nClusters=unlist(lapply(temp,'[[',1)),
    KLs=lapply(temp,'[[',2),
    cluster0=lapply(temp,'[[',3)
  )
}

KL.pick <- function(nIter, bAMatrix,Q,regQ, do.plot, seed){
  nSub <- nrow(bAMatrix)
  KL <- numeric(nSub)
  nCluster <- 1
  clusterVec <- rep(1,nSub)
  dK <- dK.func0(nCluster, clusterVec, Q, bAMatrix,regQ)
  KL[1] <- KLK.func(dK,clusterVec,bAMatrix,Q)

  clusterMatrix <- matrix(NA, nrow= nSub, ncol= nSub)
  clusterMatrix[,1] <- clusterVec

  for(nCluster in seq(2,nSub-1)){
    temp <- cluster.iter(clusterVec,nIter,nCluster,Q,bAMatrix,regQ, seed,do.plot)
    dK <- temp[[2]][,,nIter]
    clusterVec <- temp[[1]][,nIter]
    KL[nCluster] <- KLK.func(dK,clusterVec,bAMatrix,Q)
    clusterMatrix[,nCluster] <- clusterVec

  }
  list(KL=KL, clusterMatrix=clusterMatrix)
}

