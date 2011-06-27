setMethod("exprs", signature = "qPCRBatch", definition = 
    function (object) assayDataElement(object, "exprs")
)

setReplaceMethod("exprs", signature = "qPCRBatch", definition = 
    function (object, value) assayDataElementReplace(object, "exprs", value)
)

setGeneric("exprs.well.order",
    function(object)
    standardGeneric("exprs.well.order")
)

setGeneric("exprs.well.order<-",
    function(object, ..., value)
    standardGeneric("exprs.well.order<-")
)

setMethod("exprs.well.order", signature = "qPCRBatch", definition = 
    function (object) assayDataElement(object, "exprs.well.order")
)

setReplaceMethod("exprs.well.order", signature = "qPCRBatch", definition = 
    function (object, value) assayDataElementReplace(object, "exprs.well.order", value)
)
