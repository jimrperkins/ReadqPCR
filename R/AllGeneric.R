if(!isGeneric("effs")){ 
  setGeneric("effs", function(object) standardGeneric("effs"))
}

if(!isGeneric("effs<-")){ 
  setGeneric("effs<-", function(object, value) standardGeneric("effs<-"))
}

if(!isGeneric("se.effs")){ 
  setGeneric("se.effs", function(object) standardGeneric("se.effs"))
}

if(!isGeneric("se.effs<-")){ 
  setGeneric("se.effs<-", function(object, value) standardGeneric("se.effs<-"))
}

