#' Benchmark clustering output
#'
#' Evaluate clustering result by calculating Rand Index and adjusted Rand
#' Index, benchmarking against 6 traditional clustering techniques:
#' (1) HC_dist: hierarchical clustering
#' based on an inte- grated periodogram-based method as dissimilarity
#' measure (Montero and Vilar, 2014)
#' (2) HC_pred: hierarchical clustering based on a prediction
#' density-based method as dissimilarity measure (Montero and Vilar, 2014)
#' (3) BHC: Bayesian model-based hierarchical clustering with accounting
#' for uncertainty using the Dirichlet process (Savage et al., 2009)
#' (4) Mclust: finite Gaussian mixture model under Bayesian framework
#' estimated by Estimation-Maximisation (Scrucca et al., 2016)
#' (5) VC: clustering based on Bayesian mixtures of linear mixed models
#' estimated via variational inference (Tan and Nott, 2014).
#' (Omitted) KML: K-means for longitudinal data (Genolini and Falissard, 2011).
#' The relevant code is commented out but executable under local environment
#'
#' @param dat Data frame of dataset, or list of a few datasets to benchmark
#' @inheritParams clustBoot
#' @inheritParams plotPairs
#' @param y_var Character of column name in \code{dat} as observation
#' @param arr_cluster Array of predicted clusters by projection clustering method,
#' passed from function \code{\link{postPairs}}
#' @param seed Random seed to pass to benchmark methods
#'
#' @return List of matrices with rows indicating the benchmark methods and two columns
#' Rand Index and adjusted Rand Index, corresponding to the list of \code{dat}
#' @import fossil
#' @import BHC
#' @import TSclust
#' @import kml
#' @import mclust
#' @importFrom tidyr pivot_wider
#' @importFrom stats as.hclust cutree hclust
#' @export
#'
#' @examples
#' benchClust(DATASET,'ID','t','Record','Group',out_pc$arr_cluster)
#'
benchClust <- function(dat,id_var, t_var, y_var, g_var, arr_cluster, seed=1){
  if(is.data.frame(dat)) ls_dat <- list(dat)
  else ls_dat <- dat
  ls_tb <- list()

  i <- 0
  for(i in seq_along(ls_dat)){
    dat <- ls_dat[[i]]

    dat1 <- dat %>%
      select(all_of(c(id_var, t_var, y_var))) %>%
      pivot_wider(names_from = id_var, values_from= y_var) %>%
      select(-!!t_var) %>%
      as.data.frame
    dat %>% select(all_of(c(g_var, id_var))) %>% unique %>%
      select(!!g_var) %>% unlist %>% as.factor %>% as.numeric -> true_cluster
    dat2 <- t(dat1)

    nCluster <- length(unique(true_cluster))

    set.seed(seed)

    # cat('\nKML\n')
    # cld1 <<- clusterLongData(dat2)
    # kml(cld1,nCluster,10)
    # kml.c <- as.numeric(getClusters(cld1,nCluster,1))
    # fn <- "cld1.RData"
    # if (file.exists(fn))  file.remove(fn)

    cat('\nHC_dist\n')
    IP.dis <- diss(dat1, "INT.PER")#"DWT"
    IP.hclus <- cutree(hclust(IP.dis), k = nCluster)

    cat('\nHC_pred\n')
    diffs <- rep(1, ncol(dat1))
    logs <- rep(F, ncol(dat1))
    hc.hclus <- NA
    hc.hclus <- tryCatch({
      dpred <- diss(dat1, "PRED", h = nCluster, B = 1200, logarithms = logs, differences = diffs)
      cutree( hclust(dpred$dist), k = nCluster)
    },
    error=function(cond){
      message('No values returned')
      message(cond)
    })

    cat('\nBHC\n')
    if(length(unique(c(dat2)))>10){
      hc2 <- bhc(dat2, rownames(dat2), 0, seq(ncol(dat2)), "time-course",
                 numReps=1, noiseMode=0, numThreads=1, verbose=TRUE) # account for correlation in time series (must be continuous): squared exponential covariance
    }
    else hc2 <- bhc(dat2, rownames(dat2), 0, seq(ncol(dat2)), "multinomial", # used for discrete with 2+ categories data
                    numReps=1, noiseMode=0, numThreads=1, verbose=TRUE)
    bhc.hc <- cutree(as.hclust(hc2), k=nCluster)

    cat('\nMclust\n')
    BIC <- mclustBIC(dat2)
    mod1 <- Mclust(dat2, x = BIC, G = nCluster) # chosen cluster number always the minimum of G range
    mclust.c <- mod1$classification

    # listing benchmarking results
    ls_benchmark <- list(IP.hclus=IP.hclus, hc.hclus=hc.hclus, bhc.hc=bhc.hc, mclust.c=mclust.c)

    rd_ben <- sapply(list(rand.index,adj.rand.index), function(rand.fcn){
      benchmark.rand(true_cluster, ls_benchmark, rand.fcn)
    })

    cat('\nPC\n')
    rd_pc <- sapply(list(rand.index,adj.rand.index), function(rand.fcn){
      pc.rand(true_cluster, arr_cluster, rand.fcn)
    })

    tb_rand <- array(rbind(rd_ben,rd_pc), dim = c(8, 2), dimnames = list(
      c('HC_dist', 'HC_pred', 'BHC', 'Mclust',paste0('PC', seq(4))),
      c('rand','adj_rand')))
    ls_tb[[i]] <- tb_rand
  }

  ls_tb
}


benchmark.rand <- function(true_cluster, ls_benchmark, rand.fcn){
  ret <- rep(NA,length(ls_benchmark))
  for(i in seq_along(ls_benchmark)){
    if(!is.null(ls_benchmark[[i]])&length(ls_benchmark[[i]])) ret[i] <- rand.fcn(true_cluster,ls_benchmark[[i]])
  }
  ret
}

pc.rand <- function(true_cluster, arr_cluster, rand.fcn){
  apply(arr_cluster, 2, function(y)  {
    ret <- rand.fcn(true_cluster, y)
    ret
  })
}
