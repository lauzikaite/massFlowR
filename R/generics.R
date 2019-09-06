# massFlowTemplate --------------------------------------------------------
setGeneric("filepath", function(object) standardGeneric("filepath"))
setGeneric("peaksVALIDATED", function(object) standardGeneric("peaksVALIDATED"))
setGeneric("alignPEAKS", function(object, ...) standardGeneric("alignPEAKS"))
setGeneric("checkNEXT", function(object) standardGeneric("checkNEXT"))
setGeneric("validPEAKS", function(object, ...) standardGeneric("validPEAKS"))
setGeneric("fillPEAKS", function(object, ...) standardGeneric("fillPEAKS"))
setGeneric("adjustBATCH", function(object, ...) standardGeneric("adjustBATCH"))

# massFlowAnno ------------------------------------------------------------
setGeneric("annotateDS", function(object, db_file, out_dir, ...) standardGeneric("annotateDS"))
setGeneric("findANNOchemid", function(object, chemid, cutoff = 0, ...) standardGeneric("findANNOchemid"))
setGeneric("findANNOpcs", function(object, pcs, cutoff = 0, ...) standardGeneric("findANNOpcs"))
setGeneric("checkANNOTATION", function(object, chemid, pcs, out_dir = NULL, ...) standardGeneric("checkANNOTATION"))
setGeneric("comparePCS", function(object, pcs, out_dir = NULL,...) standardGeneric("comparePCS"))
setGeneric("checkADDUCTS", function(object, chemid, adducts, pcs, out_dir = NULL,...) standardGeneric("checkADDUCTS"))


