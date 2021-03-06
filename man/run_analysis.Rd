% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/driver.R
\name{run_analysis}
\alias{run_analysis}
\title{Run differential splicing anaysis.}
\usage{
run_analysis(dataDir, cases, MAX_ITER = 8, BIAS = F,
  FPKM_filter = 0.1, MCMC = FALSE, verbose = FALSE)
}
\arguments{
\item{dataDir}{A string for the path of the input data dir.}

\item{cases}{A vector of strings indicating which samples are treatment samples.}

\item{MAX_ITER}{A number for the max iteration num for the EM algorithm. Default is 8.}

\item{BIAS}{A boolean for doing bias correction or not, default is False}

\item{FPKM_filter}{A number for fpkm filtering. Default is 1.0.}

\item{MCMC}{A boolean for using MCMC instead of EM algorithm for optimization. Note that this will be very slow. Default is False.}

\item{verbose}{A boolean for more verbose messages. Default is False.}
}
\value{
A list of genes. Each gene has multiple data.frame
}
\description{
Run differential splicing anaysis.
}
\examples{
run_analysis("my_path_to_data", c("caseSample1, caseSample2, caseSample3"), MAX_ITER = 5, BIAS = F, FPKM_filter = 1.0, MCMC = FALSE, verbose = FALSE)
}
