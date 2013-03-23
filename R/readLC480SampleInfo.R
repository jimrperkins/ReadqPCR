###############################################################################
## Read Sample Information File of LC 480 into AnnotatedDataFrame
###############################################################################

read.LC480SampleInfo <- function(file, removeEmptyCols = TRUE,
                                 header = TRUE, sep = "\t", quote = "\"", 
                                 dec = ".",fill = TRUE, comment.char = "",
                                 skip = 0){
    sampleInfo <- read.table(file = file, header = header, sep = sep, 
                             quote = quote, dec = dec, fill = fill, 
                             comment.char = comment.char, skip = skip, 
                             stringsAsFactors = FALSE, check.names = FALSE)
    if(!("General:Pos" %in% names(sampleInfo)))
      stop("No 'General:Pos' (i.e. sample position) available in sample information file!")
    ## remove empty columns
    if(removeEmptyCols){
        sampleInfo <- sampleInfo[,colSums(is.na(sampleInfo)) != nrow(sampleInfo)]
    }
    ## possible column names in exported sample information file
    origNames <- c("General:Pos", "General:Sample Name", "General:Repl. Of", 
                   "General:Filt. Comb.", "General:Target Name", 
                   "General:Color", "General:Subsets", "General:Notes", 
                   "General:Sample ID", "General:Sample Prep Notes",
                   "Sample Preferences:Width", "Sample Preferences:Line Style",
                   "Sample Preferences:Point Style", "Sample Preferences:Color",
                   "Color Comp:Dominant Channel",
                   "Endpt. Geno:EndPt Sample Type", 
                   "Endpt. Geno:EndPt Genotype",
                   "Abs Quant:Sample Type", "Abs Quant:Concentration", 
                   "Abs Quant: Cp Low", "Abs Quant:Cp High", 
                   "Melt Geno:Sample Type", "Melt Geno:Genotype",
                   "Rel Quant:Target Type", 
                   "Rel Quant:Combined Sample/Target type",
                   "Rel Quant:Efficiency", 
                   "Gene Scanning:Scanning Sample Type",
                   "Gene Scanning:Scanning Genotype")
    ## new column names
    newNames <- c("Sample position", "Sample name", "Replicate of", 
                  "Filter combination", "Target name", "Color", 
                  "Subsets", "Notes", "Sample ID", "Sample Prep Notes",
                  "Sample Pref width", "Sample Pref line style", 
                  "Sample Pref point style", "Sample Pref color",
                  "Dominant channel",
                  "Endpoint sample type", "Endpoint genotype",
                  "Quantification sample type", "Concentration", 
                  "Cp Low", "Cp High", 
                  "Melt Geno sample type", "Melt Geno genotype",
                  "Target type", "Combined sample and target type", 
                  "Efficiency",                   
                  "Scanning sample type", "Scanning genotype")
    ## change column names
    for(i in 1:ncol(sampleInfo)){
        if(names(sampleInfo)[i] %in% origNames)
            names(sampleInfo)[i] <- newNames[which(origNames == names(sampleInfo)[i])]
    }
    sampleInfo[,"Sample position"] <- factor(sampleInfo[,"Sample position"], 
                                             levels = sampleInfo[,"Sample position"])
    metaData <- data.frame(labelDescription = names(sampleInfo),
                           row.names = names(sampleInfo), check.names = FALSE,
                           stringsAsFactors = FALSE)
    AnnotatedDataFrame(data = sampleInfo, varMetadata = metaData)    
}
