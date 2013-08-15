#' A deterministic statistical machine to create prognostic genomic signatures for personalized medicine. 
#'
#' The primary function used by this package is runMachine. You pass the machine: (1) a set
#' of gene expression data sets for biomarker discovery, (2) a set of gene expression data sets
#' for biomarker validation, and (3) optionally a set of gene expression data sets
#' for predictive validation. The output is a directory with a fully reproducible analysis and R
#' function for evaluating new samples. 
#' 
#' @name rdsmGeneSig
#' @author Jeff Leek <jtleek at gmail.com>
#' @references http://simplystatistics.org/2012/08/27/a-deterministic-statistical-machine/
#' 
#' @keywords dsm, personalized medicine, gene expression
#' 
#' @docType package
NULL
