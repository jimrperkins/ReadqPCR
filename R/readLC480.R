####################################################################
## Read Experiment Text File of LC 480 into CyclesSet
####################################################################

## read data from Lightcycler 480 (Roche)
read.LC480 <- function(file, colNames = c("Sample position", "Sample name",
                                          "Program number", "Segment number", 
                                          "Cycle number", "Acquisition time", 
                                          "Acquisition temperature", 
                                          "Fluorescence data"), 
                       cycleThreshold = 45, fileType = "txt", skip = 1,
                       header = TRUE, sep = "\t", quote = "\"", dec = ".",
                       fill = TRUE, comment.char = ""){
    ## txt or xml file
    if(fileType == "txt")
      res <- .readLC480txt(file = file, colNames = colNames, 
                          cycleThreshold = cycleThreshold, skip = skip,
                          header = header, sep = sep, quote = quote, dec = dec,
                          fill = fill, comment.char = comment.char)
    else
      res <- .readLC480xml(file = file, colNames = colNames, dec = dec,
                          cycleThreshold = cycleThreshold)
    res
}

## not exported function for reading xml files from LC 480
.readLC480xml <- function(file, colNames, dec, cycleThreshold){
    stop("read function for xml-files from LC 480 not yet implemented!")
}

## not exported function for reading txt files from LC 480
.readLC480txt <- function(file, colNames, cycleThreshold, header, sep, quote,
                         dec, fill, comment.char, skip){
    allData <- read.table(file = file, header = header, sep = sep, quote = quote, 
                          dec = dec, fill = fill, comment.char = comment.char, 
                          skip = skip, stringsAsFactors = FALSE)
    
    ## remove empty columns
    allData <- allData[,colSums(is.na(allData)) != nrow(allData)]
  
    ## rename of columns
    if(ncol(allData) == length(colNames))
        colnames(allData) <- colNames
    else
        stop("Number of non-empty columns not equal to length of 'colNames'!")
    
    ## split data by Sample position
    fac <- factor(allData[,"Sample position"], 
                  levels = unique(allData[,"Sample position"]))
    allData.split <- split(x = allData, f = fac)
    
    ## extract Fluorescence data
    fluoData <- sapply(allData.split, function(x) x[seq_len(cycleThreshold),"Fluorescence data"])
    rownames(fluoData) <- seq_len(cycleThreshold)
    
    ## feature Data (here: cycles)
    acTime <- allData.split[[1]][seq_len(cycleThreshold),"Acquisition time"]
    acTemp <- allData.split[[1]][seq_len(cycleThreshold),"Acquisition temperature"]
    fData <- data.frame("Cycle number" = seq_len(cycleThreshold),
                        "Acquisition time" = acTime,
                        "Acquisition temperature" = acTemp,
                        row.names = seq_len(cycleThreshold), 
                        check.names = FALSE,
                        stringsAsFactors = FALSE)
    fMetaData <- data.frame(labelDescription = c("Cycle number", 
                                                 "Acquisition time",
                                                 "Acquisition temperature"),
                            row.names = names(fData), check.names = FALSE,
                            stringsAsFactors = FALSE)
    featureData <- AnnotatedDataFrame(data = fData, varMetadata = fMetaData)
    
    ## phenotypic data
    samPos <- sapply(allData.split, function(x) x[1,"Sample position"])
    samPos <- factor(samPos, levels = unique(samPos))
    samNam <- sapply(allData.split, function(x) x[1,"Sample name"])
    prgNum <- sapply(allData.split, function(x) x[1,"Program number"])
    segNum <- sapply(allData.split, function(x) x[1,"Segment number"])
  
    pData <- data.frame("Sample position" = samPos,
                        "Sample name" = samNam,
                        "Program number" = prgNum,
                        "Segment number" = segNum,
                        row.names = samPos, check.names = FALSE,
                        stringsAsFactors = FALSE)
    pMetaData <- data.frame(labelDescription = c("Sample position", 
                                                 "Sample name",
                                                 "Program number",
                                                 "Segment number"),
                            row.names = names(pData), check.names = FALSE,
                            stringsAsFactors = FALSE)
    phenoData <- AnnotatedDataFrame(data = pData, varMetadata = pMetaData)
    
    new("CyclesSet", exprs = fluoData, featureData = featureData, 
        phenoData = phenoData)
}
