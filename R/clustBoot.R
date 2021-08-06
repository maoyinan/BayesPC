#' Choose cluster number by bootstrap method
#'
#' With predicted projections from \code{pcFit}, take bootstrap samples for
#' \code{nBoot} times, at each time two random samples are drawn and applied in
#' \code{kmeans} clustering with both input sets to create 4 sets of clustering results.
#' Empirical clustering distance is hence computed to enable choosing minimal cluster number
#' with cutoff being half of the largest clustering distance.
#'
#' @inheritParams modelStan
#' @param t_var Column name in \code{dat} representing the scaled time points ranging in \eqn{[0,1]}
#' @param mat_fitted Prediction replicates of Linear Mixed Model output. Multiple columns
#' corresponds to sets of random effects of projection defined in \code{pcFit}
#' @param nB Number of bootstrap samples
#' @param nCl Maximum number of clusters to consider
#' @param seed Random seed to pass to bootstrap
#'
#' @return list(sBtb, cutoff, nClusters)
#' \item{sBtb}{Matrix of clustering distance, with \code{nCl} rows and \code{ncol(mat_fitted} columns
#' corresponding to projections in \code{mat_fitted}}
#' \item{cutoff}{Vector of cutoff values: \eqn{cutoff_j=1/2*max(sBtb[,j])}}
#' \item{nClusters}{Vector of cluster numbers chosen by the minimum number with
#' clustering distance no larger than the cutoff}
#'
#' @family Cluster number choices
#' @importFrom graphics text
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' data(df_of_draws)
#' ls_par <- postMean(df_of_draws, paste0("Z", 1:10), "ID", DATASET)
#' ls_idxA <- list(
#'   seq(10),
#'   1:4,
#'   5:7,
#'   8:10
#' )
#' mat_fitted <- pcFit("ID", ls_par, DATASET, ls_idxA)
#' out_boot <- clustBoot(paste0("Z", 1:10), "ID","t", mat_fitted,DATASET, 10, 5, 1)

clustBoot <- function(x_var, id_var, t_var, mat_fitted, dat, nB, nCl, seed){
  set.seed(seed)
  nSet <- ncol(mat_fitted)
  sBtb <- matrix(NA,nrow=nCl-1, ncol=nSet)
  nClusters = cutoff <- numeric(nSet)
  y <- 0

  for(j in seq(nSet)){
    xx <- dat %>%
      dplyr::select(!!id_var,!!t_var) %>%
      dplyr::bind_cols(data.frame(y=mat_fitted[,j])) %>%
      tidyr::pivot_wider(names_from = !!t_var, values_from=y) %>%
      dplyr::select(-!!id_var) %>%
      as.matrix
    nSub <- nrow(xx)
    sB <- sapply(seq(2,nCl), function(cl){
      cat('\nSet',j,'Bootstrap with',cl,'clusters')
      mean(sapply(seq(nB),function(b){
        idx1 <- sample(nSub,nSub,replace=T)
        idx2 <- sample(nSub,nSub,replace=T)
        out1 <- stats::kmeans(xx[idx1,], cl)
        out2 <- stats::kmeans(xx[idx2,], cl)
        fitted1 <- stats::fitted(out1, method='classes')
        fitted12 <- my.dist(xx[idx2,], out1$centers)
        fitted2 <- stats::fitted(out2, method='classes')
        fitted21 <- my.dist(xx[idx1,], out2$centers)

        # empirical clustering distance
        idx_pair <- matrix(c(rep(1:nSub,nSub),rep(1:nSub,each=nSub)),nrow=2, ncol=nSub^2, byrow = T)
        d <- mean(apply(idx_pair,2,function(x){
          (fitted1[x[1]]==fitted12[x[2]] & fitted21[x[1]]!=fitted2[x[2]] )+
            ( fitted1[x[1]]!=fitted12[x[2]] & fitted21[x[1]]==fitted2[x[2]])
        }))
      }))
    })
    sBtb[,j] <- sB
    cutoff[j] <- max(sB)/2
    idx <- setdiff(which(sB<=cutoff[j]), 1)[1]
    nClusters[j] <- idx+1
  }
  ls_sB <- list(sBtb=sBtb, cutoff=cutoff, nClusters=nClusters)
  ls_sB
}

# helper function
my.dist <- function(x, centers) {
  # compute squared euclidean distance from each sample to each cluster center
  tmp <- sapply(seq_len(nrow(x)),
                function(i) apply(centers, 1,
                                  function(v) sum((x[i, ]-v)^2)))
  if(!is.matrix(tmp)) tmp <- t(tmp)
  max.col(-t(tmp))  # find index of min distance
}

