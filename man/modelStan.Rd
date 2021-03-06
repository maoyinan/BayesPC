% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelStan.R
\name{modelStan}
\alias{modelStan}
\title{Bayesian mixed linear regression modelling}
\usage{
modelStan(y_var, x_var, id_var = "ID", dat, seed = 1, ...)
}
\arguments{
\item{y_var}{Character of response variable}

\item{x_var}{Character vector of random effect variables}

\item{id_var}{Character of id variable}

\item{dat}{Longitudinal data input}

\item{seed}{Random seed to pass to \link[rstan]{sampling}}

\item{...}{Additional arguments passed to \link[rstan]{sampling}}
}
\value{
dataframe of posterior output from Bayesian mixed linear regression
}
\description{
To fit longitudinal dataset using linear mixed model by calling stan package.
}
\examples{
df_of_draws <- modelStan("Record", paste0("Z", 1:10), "ID", DATASET, seed=1, iter = 2, chain = 1)
}
\seealso{
Other BayesPC main functions: 
\code{\link{clustKL}()},
\code{\link{pcFit}()},
\code{\link{postMean}()}
}
\concept{BayesPC main functions}
