setGeneric("filepath", function(object) standardGeneric("filepath"))
setGeneric("peaksVALIDATED", function(object) standardGeneric("peaksVALIDATED"))
setGeneric("alignPEAKS", function(object, out_dir, ncores = 2, ...) standardGeneric("alignPEAKS"))
setGeneric("checkNEXT", function(object) standardGeneric("checkNEXT"))
setGeneric("validPEAKS", function(object, out_dir,  ncores = 2, ...) standardGeneric("validPEAKS"))
setGeneric("fillPEAKS", function(object, fill_value = "into", out_dir, ncores = 2, ...) standardGeneric("fillPEAKS"))
setGeneric("adjustBATCH", function(object, out_dir, batch_next_file, batch_end_roi, batch_start_roi, ...) standardGeneric("adjustBATCH"))


