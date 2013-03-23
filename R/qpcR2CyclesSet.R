###############################################################################
## Transform data sets of package qpcR to CyclesSet
###############################################################################

qpcR2CyclesSet <- function(x, cyc = 1, cycleThreshold){
  x <- x[order(x[,cyc]),]
  if(missing(cycleThreshold)){
      fluoData <- x[,-cyc]
      rownames(fluoData) <- seq_len(nrow(fluoData))
  }else{
      fluoData <- x[seq_len(cycleThreshold),-cyc]
      rownames(fluoData) <- seq_len(cycleThreshold)
  }
  fData <- data.frame("Cycle number" = seq_len(nrow(fluoData)),
                      row.names = seq_len(nrow(fluoData)), 
                      check.names = FALSE,
                      stringsAsFactors = FALSE)
  fMetaData <- data.frame(labelDescription = "Cycle number",
                          row.names = names(fData), check.names = FALSE,
                          stringsAsFactors = FALSE)
  featureData <- new("AnnotatedDataFrame", data = fData, varMetadata = fMetaData)
  
  samNam <- colnames(fluoData)
  
  pData <- data.frame("Sample name" = samNam,
                      row.names = samNam, check.names = FALSE,
                      stringsAsFactors = FALSE)
  pMetaData <- data.frame(labelDescription = "Sample name",
                          row.names = names(pData), check.names = FALSE,
                          stringsAsFactors = FALSE)
  phenoData <- new("AnnotatedDataFrame", data = pData, varMetadata = pMetaData)
  
  new("CyclesSet", exprs = fluoData, featureData = featureData, 
      phenoData = phenoData)  
}

