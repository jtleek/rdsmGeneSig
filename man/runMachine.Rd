\name{runMachine}
\alias{runMachine}
\title{Run the deterministic statistical machine}
\usage{
  runMachine(discCEL, discMETA, validCEL = NULL,
    validMETA = NULL,
    reportTitle = "Gene Signature Development",
    dirName = "geneSigOutput", openReport = TRUE,
    analysisSeed = 32313)
}
\arguments{
  \item{discCEL}{A character vector where entries are
  directories with the CEL files for development prognostic
  data sets}

  \item{discMETA}{A character vector with entries being csv
  files for meta data for development prognostic data sets}

  \item{validCEL}{(optional) AA character vector where
  entries are directories with the CEL files for validation
  data sets}

  \item{validMETA}{(optional) A character vector with
  entries being csv files for validation data sets}

  \item{predCEL}{(optional) AA character vector where
  entries are directories with the CEL files for predictive
  data sets}

  \item{predMETA}{(optional) A character vector with
  entries being csv files for predictive data sets}

  \item{dirName}{The base name of the directory to be
  created}

  \item{openReport}{If TRUE, the HTML report is opened in
  the browser.}

  \item{seed}{Set the seed for the analysis.}
}
\value{
  The directory with all of the output of the DSM.
}
\description{
  This function runs the deterministic statistical machine
  on a set of discovery and validation data sets and
  optionally on a clinical data set for prediction
  validation.
}

