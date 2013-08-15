#' Run the deterministic statistical machine
#' 
#' This function runs the deterministic statistical machine on a set of discovery
#' and validation data sets and optionally on a clinical data set for prediction
#' validation. 
#' 
#' @param discCEL A character vector where entries are directories with the CEL files for development prognostic data sets
#' @param discMETA A character vector with entries being csv files for meta data for development prognostic data sets
#' @param validCEL (optional) AA character vector where entries are directories with the CEL files for validation data sets
#' @param validMETA (optional) A character vector with entries being csv files for validation data sets
#' @param predCEL (optional) AA character vector where entries are directories with the CEL files for predictive data sets
#' @param predMETA (optional) A character vector with entries being csv files for predictive data sets
#' @param dirName The base name of the directory to be created 
#' @param openReport If TRUE, the HTML report is opened in the browser.
#' @param seed Set the seed for the analysis. 
#' 
#' @export
#' 
#' @return The directory with all of the output of the DSM. 

runMachine <- function(discCEL, discMETA, 
                       validCEL=NULL,validMETA=NULL, 
                       reportTitle="Gene Signature Development",
                       dirName="geneSigOutput", 
                       openReport=TRUE,analysisSeed=32313){
  
  require(knitrBootstrap)
  
  # Create the directory
  dir.create(dirName)
  
  setwd(dirName)
  # Set up the file names
  mdFileName <- paste(dirName,".md",sep="")
  htmlFileName <-  paste(dirName,".html",sep="")
  
  # Knit the file
  workDir <- tools:::file_path_as_absolute(paste0("../",dirName))
    
  machineRmd <- system.file(package="rdsmGeneSig","geneSig.Rmd")
 
  knit2html(input=machineRmd,output=mdFileName)
  bootstrap_HTML(input=htmlFileName,output =paste0("bootstrap",htmlFileName))
  if(openReport){browseURL(paste0("bootstrap",htmlFileName))}
  setwd("../")
}