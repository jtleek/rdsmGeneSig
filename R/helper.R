makeBatch <- function(object){
  
  nSamples <- dim(pData(object))[1]
  
  tmp1 <- sum(grepl("batch",tolower(names(pData(object)))))
  if(sum(tmp1)==1){
    pdBatch <- pData(object)[,grepl("batch",tolower(names(pData(object))))]
  }else{
    pdBatch <- rep(NA,nSamples)
  }
  
  tmp2 <- sum(grepl("rin",tolower(names(pData(object)))))
  if(sum(tmp1)==1){
    rin <- pData(object)[,grepl("rin",tolower(names(pData(object))))]
  }else{
    rin <- rep(NA,nSamples)
  }
  
  dates=vector("character",nSamples)
  files <- as.character(pData(object)$fullFilePath)
  for(i in seq(along=dates)){
    tmp=affyio::read.celfile.header(files[i],info="full")$ScanDate
    if(length(tmp) > 0){
      dates[i]=strsplit(tmp," ")[[1]][1]
    }
  }
  scanBatch <- make.consecutive.int(as.factor(dates))
  if(allequal(scanBatch)){scanBatch <- rep(NA,nSamples)}
  if(allequal(pdBatch)){pdBatch <- rep(NA,nSamples)}
  return(list(scanBatch=scanBatch,pdBatch=pdBatch,rin=rin))
}

unzipCels <- function(dir){
  fileList <- list.files(dir,full.names=TRUE)
  fileList <- fileList[grepl(".gz",fileList)]
  for(i in 1:length(fileList)){
    gunzip(fileList[i])
  }
}

make.consecutive.int <- function(y) {
  oldWarn = getOption("warn")
  ## Turn off warnings.
  options(warn = -1)
  
  if(is.null(y)) {return(NULL)}
  
  if(!is.vector(y))
    y = as.vector(as.character(y))
  
  out <- as.integer(as.factor(as.character(y)))
  
  options(warn = oldWarn)
  
  return(out)
}


getPheno <- function(object){
  tmp1 <- match("outcome",tolower(names(pData(object))))
  if(!is.na(tmp1) & length(tmp1)==1){
    outcome <- as.numeric(pData(object)[,tmp1])
  }else{
    print("You must have a unique variable called outcome in your metadata.")
  }
 

  tmp1 <- match("survivaltime",tolower(names(pData(object))))
  if(!is.na(tmp1) & length(tmp1)==1){
    survivalTime <- as.numeric(pData(object)[,tmp1])
  }else{
    print("You must have a unique variable called survival time in your metadata.")
  }
  
  outcome <- Surv(survivalTime,outcome)
  
  return(outcome)
}

allequal <- function(x){
  if(length(unique(x))==1){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

batchOutText <- function(vv,oo,i){
  if(class(oo)=="Surv"){
    pp <- anova(coxph(oo ~ vv))[2,4]
  }else{
    pp <- anova(lm(oo ~ as.factor(vv)))[1,5]
  }
  
  if(i == 1){
    tt <- paste("We determined that there was an [RNA quality (RIN)](http://en.wikipedia.org/wiki/RNA_integrity_number) variable in your metadata. The P-value for the association between RIN and your outcome was",round(pp,3),".")
    if(pp > 0.05){
      tt<- paste(tt,"This suggests there is not a strong association between outcome and RIN.")
    }else{
      tt <- paste(tt,"This suggests there is a strong association between your outcome and RIN. This is a potential problem with your experimental design and predictions should be used cautiously.Please consult a statistician about potential study design issues.")
    }
  }
  if(i==2){
    tt <- paste("We determined that there was a batch  variable in your metadata. The P-value for the association between your batch variable and your outcome was",round(pp,3),".")
    if(pp > 0.05){
      tt<- paste(tt,"This suggests there is not a strong association between outcome and your batch variable.")
    }else{
      tt <- paste(tt,"This suggests there is a strong association between your outcome and your batch variable. This is a potential problem with your experimental design and predictions should be used cautiously.Please consult a statistician about potential study design issues.")
    }
  }
  if(i==3){
    tt <- paste("We determined calculated a batch variable based on the scan dates of your microarrays. The P-value for the association between scan batch and your outcome was",round(pp,3),".")
    if(pp > 0.05){
      tt<- paste(tt,"This suggests there is not a strong association between outcome and scan batch.")
    }else{
      tt <- paste(tt,"This suggests there is a strong association between your outcome and scan batch. This is a potential problem with your experimental design and predictions should be used cautiously.Please consult a statistician about potential study design issues.")
    }
  }
  return(tt)
}


makeFormula <- function(pairMat){
  tmp <- paste("phenoVar~")
  tmp2 <- paste(paste0("'",colnames(pairMat),"'"),collapse="+")
  return(formula(paste(tmp,tmp2)))
}

quantileSequence <- function(x,n){
  return(seq(quantile(x,0.2),quantile(x,0.8),length=n))
}

discMergeText <-"The annotation files are merged using an [outer join](http://stackoverflow.com/questions/8091303/merge-multiple-data-frames-in-a-list-simultaneously) and the expression sets are merged using the protocal from [inSilicoMerging](http://www.bioconductor.org/packages/2.12/bioc/html/inSilicoMerging.html)." 

discMergeText2 <-"There is no merging necessary."

