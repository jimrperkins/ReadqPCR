setMethod("exprs", signature = "CyclesSet", definition = 
            function (object) assayDataElement(object, "exprs")
)

setReplaceMethod("exprs", signature = "CyclesSet", definition = 
                   function (object, value) assayDataElementReplace(object, "exprs", value)
)
