% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evluate.R
\name{evaluate_gene_levels}
\alias{evaluate_gene_levels}
\title{Compare Raspberry gene-level results against other tools.}
\usage{
evaluate_gene_levels(raspberry_result, truth_file, ...)
}
\arguments{
\item{raspberry_result}{the results from run_analysis.}

\item{truth_file}{A file contains differential spliced gene id. One at each line.}

\item{...}{Other programs output as a key value pari, separated by comma. The keys are program names, which have to be from the following (CUffdiff2, DEXSeq, DSGseq. currently supported). The values are file path to the program  output.}
}
\value{
A path the ROC curve pdf.
}
\description{
Compare Raspberry gene-level results against other tools.
}
\examples{
evaluate_gene_levels(raspberry$result, truth_file = "inst/extdata/trueASgenes.txt",
                    Cuffdiff2 = "inst/extdata/splicing.diff",
                    DEXSeq = "inst/extdata/dexseq.txt",
                    DSGseq = "inst/extdata/DSGresult.txt")
}
