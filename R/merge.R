###############################################################################
## methods for merge
###############################################################################

## make function merge generic
if(!isGeneric("merge")){ 
  setGeneric("merge", function(x, y, ...) standardGeneric("merge"))
}

## make function merge generic
setMethod("merge", signature(x = "eSet", 
                             y = "AnnotatedDataFrame"),
          function(x, y, eSet.slot = "phenoData", 
                   by = intersect(names(pData(x)), names(pData(y))),
                   by.x = by, by.y = by, all = FALSE, all.x = all, all.y = all,
                   sort = FALSE, suffixes = c(".x",".y"), incomparables = NULL, ...){
              pdat <- merge(pData(x), pData(y), by = by, by.x = by.x, 
                            by.y = by.y, all = all, all.x = all.x, 
                            all.y = all.y, sort = sort, 
                            suffixes = suffixes, incomparables = NULL, ...)
              if(eSet.slot == "phenoData"){
                  metdat <- rbind(varMetadata(phenoData(x)), varMetadata(y))
                  metdat <- metdat[!duplicated(metdat), ,drop = FALSE]
                  pData(x) <- pdat
                  varMetadata(phenoData(x)) <- metdat
              }
              if(eSet.slot == "featureData"){
                  metdat <- rbind(varMetadata(featureData(x)), varMetadata(y))
                  metdat <- metdat[!duplicated(metdat), ,drop = FALSE]
                  fData(x) <- pdat
                  varMetadata(featureData(x)) <- metdat
              }
              if(eSet.slot == "protocolData"){
                  metdat <- rbind(varMetadata(protocolData(x)), varMetadata(y))
                  metdat <- metdat[!duplicated(metdat), ,drop = FALSE]
                  pData(protocolData(x)) <- pdat
                  varMetadata(protocolData(x)) <- metdat
              }
              x
          })

setMethod("merge", signature(x = "AnnotatedDataFrame", 
                             y = "eSet"),
          function(x, y, eSet.slot = "phenoData", 
                   by = intersect(names(pData(x)), names(pData(y))),
                   by.x = by, by.y = by, all = FALSE, all.x = all, all.y = all,
                   sort = FALSE, suffixes = c(".x",".y"), incomparables = NULL, ...){
          merge(x=y, y=x, eSet.slot = eSet.slot, by = by, by.x = by.x, by.y = by.y,
                all = all, all.x = all.x, all.y = all.y, sort = sort, 
                suffixes = suffixes, incomparables = incomparables, ...)
          })