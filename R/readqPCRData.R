read.taqman <- function(..., filenames = character(0), phenoData = new("AnnotatedDataFrame"), notes = "", verbose = FALSE)
{
    auxnames <- unlist(list(...))
    filenames <- c(filenames, auxnames)
    checkValidTaqmanFilenames(filenames)
    pdata <- pData(phenoData) # number of files
    taqInfo <- .read.TaqBatch(filenames, verbose) # need to make this work for tech reps and multiple files
    exprs <- taqInfo$exprs
    well.order <- taqInfo$well.order
    exprs.well.order <- assayDataNew("environment", exprs.well.order = well.order)
    n <- length(colnames(exprs))
    if (dim(pdata)[1] != n) { # so if we don't have a row for each sample in the pData matrix
        warning("Incompatible phenoData object. Created a new one using sample name data derived from raw data.\n")
        samplenames <- sub("^/?([^/]*/)*", "", colnames(exprs))
        pdata <- data.frame(sample = 1:length(samplenames), row.names = samplenames)
        phenoData <- new("AnnotatedDataFrame", data = pdata,
            varMetadata = data.frame(labelDescription = "arbitrary numbering",
                row.names = "sample"))
    }
    return(new("qPCRBatch", exprs = exprs, phenoData = phenoData, exprs.well.order = well.order))
}

.read.TaqBatch <- function(filenames, verbose)
{
    totalPlateIds <- vector()
    plate.offset <- 0 # this is used to name plates when combining several files 
    fileNameCount <- 1 # used to count number of files given
    for (filename in filenames) {
        if (verbose) message("Filename[ ",filename," ]")
        raw.data <- read.delim(filename, skip = 12) # read in file, ignoring the gubbins on the first 12 lines
        if (! 1 %in% regexpr("Summary", as.character(raw.data[,1]))) stop("Problems with Taqman file, Summary info not found") #
        EndOfData <- grep("Summary", raw.data[,1])
        raw.data <- raw.data[1:EndOfData-1, ] # get rid of from where data finishes until the end of the file


        raw.data$Detector <- make.names(raw.data$Detector)
        raw.data$Sample <- make.names(raw.data$Sample)
        samples <- unique(raw.data$Sample) # individual sample names
        detectors <- unique(raw.data$Detector) # individual detector names

        if(fileNameCount > 1) { # do some checking if we are the second time around
          if(TRUE %in% (samples == colnames(totalExprs))) stop("Can't combine files, > 1 sample labels are the same between samples")
          if(TRUE %in% (raw.data$PlateID %in% totalPlateIds)) 
            stop ("Can't proceed, duplicate plate Ids in different files. All plate IDs should be unique")
        }
        original.order <- list() # initialise the list to keep track of the original order of the 
	raw.data$Ct[as.character(raw.data$Ct) %in% "Undetermined"] <- NA
        ############################################################
        # Add Plate ID information IF there were none to being with
        if ("" %in% as.character(raw.data$PlateID)) {
          wells.per.plate <- max(as.numeric(as.character(raw.data$Well)))
          number.of.plates <- length(as.character(raw.data$Well)) / wells.per.plate
          well.names <- vector(length = length(as.character(raw.data$Well)))
          plate.well <- 1
          for (plate.number in 1:number.of.plates) {
            for (well.number in 1:wells.per.plate) {
               well.names[plate.well] <- paste(plate.number + plate.offset, well.number, sep= "-")
               plate.well <- plate.well + 1
            }
          }
          plate.offset <- plate.offset + number.of.plates # add the number of plates to offset the plate number
          raw.data$PlateID <- well.names
        }
        else { ## otherwise we just paste the well number onto the Plate ID value
          totalPlateIds <- raw.data$PlateID # put them in a variable for checking for duplication
          raw.data$PlateID <- paste(raw.data$PlateID, as.character(raw.data$Well), sep= "-")
        }
        ####################################################
        allDetectors <- raw.data$Detector # save these as a variable because raw.data$Detector will get changed
        firstTimeFlag <- TRUE
        for (sample in samples) { # for each sample
            if (verbose) message("Now reading for sample:", sample, "\n")
            total.detectors <- length(allDetectors[raw.data$Sample == sample])
            individual.detectors <- length(unique(allDetectors[raw.data$Sample == sample]))
            tech.reps <- total.detectors/individual.detectors
            if ((tech.reps %% 1) != 0) { # if total number of replicates not a multiple of number of individual detectors
              warning.text <- paste("Corrupt taqman file: total number of readings for sample ", 
                sample, " not a multiple of number of individual number of detectors")
              stop(warning.text)
            }
            if (tech.reps > 1) { # Currently can't cope with technical replicates
              if(verbose) message("More than 1 technical replicate detected\n")
              staticDetector <- raw.data$Detector[raw.data$Sample == sample] # need this since we are modifying raw.data$Detector
              for(techDetect in unique(raw.data$Detector[raw.data$Sample == sample]))  {
                techDLength <- sum(staticDetector %in% techDetect)
                suffixedNames <- paste(techDetect, 1:techDLength, sep="_TechReps.")
                raw.data$Detector[raw.data$Sample == sample][raw.data$Detector[raw.data$Sample == sample] %in% techDetect] <- suffixedNames
              }
            }
            if(fileNameCount > 1) { # do some checking
              if(! FALSE %in% (sort(raw.data$Detector[raw.data$Sample == sample]) == sort(rownames(totalExprs)))) {
                if(verbose == TRUE) message("we are combining files with the same detector names\n")
              }
              else stop("Problem combining files on detector names. Make sure detector names match for all files\n")
              }#raw.data$Detector <- as.character(raw.data$Detector) # coerce to stop funny behaviour
              if(firstTimeFlag == TRUE) {
              exprs <- data.frame(unique(raw.data$Detector), row.names=1) # start the exprs data frame
              well.order <- data.frame(unique(raw.data$Detector), row.names=1)
              firstTimeFlag <- FALSE
              }
            original.order <- c(original.order,list(cbind(as.character(raw.data$Detector[raw.data$Sample == sample]),
                as.character(raw.data$Ct[raw.data$Sample == sample])))) # This bit to add the information about pipetting and order

            well.info <- data.frame(raw.data$Detector[raw.data$Sample == sample], # put Well data in a matrix
                         raw.data$PlateID[raw.data$Sample == sample],
                         row.names=1)
            Cts <- data.frame(raw.data$Detector[raw.data$Sample == sample], # put Cts values in a matrix
                         as.numeric(as.character(raw.data$Ct[raw.data$Sample == sample])),
                         row.names=1)
            exprs <- data.frame(merge(exprs, Cts, by="row.names"), row.names=1)
            if (! FALSE %in% (row.names(exprs) == row.names(Cts))) stop("BYE")

            well.order <- data.frame(merge(well.order, well.info, by="row.names"), row.names=1)
            if (verbose) message("sample ", sample, "read\n")
        }
        names(well.order) <- samples
        names(exprs) <- samples
        exprs <- as.matrix(exprs)
        well.order <- as.matrix(well.order)
        if(fileNameCount == 1) {
          totalExprs <- exprs
          totalWell.order <- well.order
        }
        else {
          totalWell.order <- cbind(totalWell.order, well.order)
          totalExprs <- cbind(totalExprs,exprs)
        }
        fileNameCount <- fileNameCount + 1
    }
    taqInfo <- list()
    taqInfo$exprs <- totalExprs
    taqInfo$well.order <- totalWell.order
    return(taqInfo)
}

checkValidTaqmanFilenames <- function (filenames)
{
    if (!is.character(filenames))
        stop(strwrap(paste("file names must be specified using a character",
            "vector, not a", sQuote(typeof(filenames)))), call. = FALSE)
    if (length(filenames) == 0)
        stop("no file names provided")
    if (any(sapply(filenames, nchar) < 1))
        stop("empty file names are not allowed")
    finfo <- file.info(filenames)
    whBad <- sapply(finfo[["isdir"]], function(x) !identical(FALSE,
        x))
    if (any(whBad)) {
        msg <- paste("the following are not valid files:\n",
            paste("  ", filenames[whBad], collapse = "\n"))
        stop(msg, call. = FALSE)
    }
    TRUE
}

read.qPCR <- function(filename = character(0), phenoData = new("AnnotatedDataFrame"), notes = "", verbose = FALSE)
{
    pdata <- pData(phenoData)
    checkValidqPCRFilename(filename)
    qPCRInfo <- .read.qPCR(filename, verbose) # need to make this work for tech reps and multiple files
    exprs <- qPCRInfo$exprs
    well.order <- qPCRInfo$well.order
    exprs.well.order <- assayDataNew("environment", exprs.well.order = exprs)
    n <- length(colnames(exprs))
    if (dim(pdata)[1] != n) { # so if we don't have a row for each sample in the pData matrix
        warning("Incompatible phenoData object. Created a new one using sample name data derived from raw data.\n")
        samplenames <- sub("^/?([^/]*/)*", "", colnames(exprs))
        pdata <- data.frame(sample = 1:length(samplenames), row.names = samplenames)
        phenoData <- new("AnnotatedDataFrame", data = pdata,
            varMetadata = data.frame(labelDescription = "arbitrary numbering",
                row.names = "sample"))
    }
    if(! is.null(qPCRInfo$well.order)) {
        return(new("qPCRBatch", exprs = exprs, phenoData = phenoData, exprs.well.order = well.order))
    }
    else {

        return(new("qPCRBatch", exprs = exprs, phenoData = phenoData))
    }
}

.read.qPCR <- function(filename, verbose)
{
    noWellData <- FALSE
    raw.data <- read.table(filename, header=TRUE)
    if(is.null(raw.data$Well) || is.null(raw.data$PlateID)) {
         noWellData <- TRUE
         if (verbose) message("No Well and/or Plate info found, skipping this part", "\n")
    }
    else {
        raw.data$PlateID <- paste(raw.data$PlateID, as.character(raw.data$Well), sep= "-")
    }
    levels(raw.data$Sample) <- make.names(levels(raw.data$Sample))
    levels(raw.data$Detector) <- make.names(levels(raw.data$Detector))
    Ct <- as.character(raw.data$Ct)
    samples <- levels(raw.data$Sample)
    detectors <- levels(raw.data$Detector)
    allDetectors <- raw.data$Detector
    firstTimeFlag <- TRUE
    for (sample in samples) { # for each sample
        if (verbose) message("Now reading for sample:", sample, "\n")
        total.detectors <- length(allDetectors[raw.data$Sample == sample])
        individual.detectors <- length(levels(allDetectors[raw.data$Sample == sample]))
        tech.reps <- total.detectors/individual.detectors
        raw.data$Detector <- as.character(raw.data$Detector)
          if ((tech.reps %% 1) != 0) { # if total number of replicates not a multiple of number of individual detectors
            warning.text = paste("File incorrect, make sure that detectors are the same for all samples")
            stop(warning.text)
        }
        if (tech.reps > 1) {
            if(verbose) message("More than 1 technical replicate detected\n")
              staticDetector <- raw.data$Detector[raw.data$Sample == sample]
              for(techDetect in unique(raw.data$Detector[raw.data$Sample == sample]))  {
                techDLength <- sum(staticDetector %in% techDetect)
                suffixedNames <- paste(techDetect, 1:techDLength, sep="_TechReps.")
                raw.data$Detector[raw.data$Sample == sample][raw.data$Detector[raw.data$Sample == sample] 
                  %in% techDetect] <- suffixedNames
            }
       }
       if(firstTimeFlag == TRUE) {
           exprs <- data.frame(unique(raw.data$Detector), row.names=1) # start the exprs data frame
           well.order <- data.frame(unique(raw.data$Detector), row.names=1)
           firstTimeFlag <- FALSE
       }
       raw.data$Detector <- as.factor(raw.data$Detector)
        if(noWellData == FALSE) {
            well.info <- data.frame(raw.data$Detector[raw.data$Sample == sample], # put well info values in a matrix
              raw.data$PlateID[raw.data$Sample == sample],
                row.names=1)
        }
        Cts <- data.frame(raw.data$Detector[raw.data$Sample == sample],
          as.numeric(as.character(raw.data$Ct[raw.data$Sample == sample])),
            row.names=1)
        exprs <- data.frame(merge(exprs, Cts, by="row.names"), row.names=1)
        if(noWellData == FALSE) {
            well.order <- data.frame(merge(well.order, well.info, by="row.names"), row.names=1)
        }
        if (verbose) message("sample ", sample, "read\n")
    }
    qPCRInfo <- list()
    names(exprs) <- samples
    qPCRInfo$exprs <- as.matrix(exprs)
    if(noWellData == FALSE) {
      colnames(well.order) <- names(exprs)
      qPCRInfo$well.order <- well.order
    }
    return(qPCRInfo)
}

checkValidqPCRFilename <- function (filename)
{
    if (!is.character(filename))
        stop(strwrap(paste("file name must be specified using a character",
            "vector, not a", sQuote(typeof(filename)))), call. = FALSE)
    if (length(filename) == 0)
        stop("no file name provided")
    if (any(sapply(filename, nchar) < 1))
        stop("empty file name not allowed")
    finfo <- file.info(filename)
    whBad <- sapply(finfo[["isdir"]], function(x) !identical(FALSE,
        x))
    if (any(whBad)) {
        msg <- paste("not valid file:\n",
            paste("  ", filename[whBad], collapse = "\n"))
        stop(msg, call. = FALSE)
    }
    TRUE
}
