#' Sample \code{\link{postPairs}} output to demonstrate use of package.
#'
#' List of 2 elements.
#' 1)  A list containing matrix of probabilities respective to random effect sets,
#' with rows and columns corresponding to subject ID. Cell values lie between 0 and 1.
#' Co-membership probability table is based on 1000 posterior draws based on output
#' from \code{\link{postPairs}} function using the simulated dataset \code{DATASET},
#' accessible by \code{data(DATASET)}.
#' 2) Array of cluster labels based on the posterior draws, with each column
#' corresponding to given random effects.
#'
#' @format list
#'
"out_pc"
