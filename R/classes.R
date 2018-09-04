setClass("massFlowDB",
         slots = c(files = "data.frame",
                   db = "data.frame"
         ))
setClass("massFlowTemplate",
         slots = c(files = "list",
                   tmp = "data.frame"))

