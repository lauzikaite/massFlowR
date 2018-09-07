setClass("massFlowDB",
         slots = c(
           db_filepath = "character",
           db = "data.frame"
         ))
setClass("massFlowTemplate",
         slots = c(files = "list",
                   tmp = "data.frame"))

