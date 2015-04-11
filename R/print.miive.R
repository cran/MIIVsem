#' @export
print.miive<- function(x, digits = 3,...){
  #options(scipen=10, digits=digits)
  cat("\n")
  cat("Parameter Estimates \n")
  cat("\n")
  
  dat   <- x$dat
  model <- x$model  #paste("t", i, )
  modeqns <- x$modeqns
  `$`(dat , "P(|Z|)") <- round(`$`(dat , "P(|Z|)"),digits)
  dat$DV[duplicated(dat$DV)] <- NA
  cf <- format(dat, digits = digits) ## use format to set other options like digits, justify , ...
  cf[is.na(dat)] <- ""
  print(cf, 
      quote = FALSE, right = FALSE, row.names = FALSE,print.gap=2, na.print = "")
  
  z_ind <- 1
  for (i in 1:length(model)){
    if (model[[i]]$NOTE != "") {
      if (z_ind ==1){
          cat("\n")
          cat("Notes:")
          cat("\n")
      }
      cat(model[[i]]$NOTE, "\n")
      z_ind = z_ind + 1
    }
  }
  
  cat("\n")
  cat("\n")
  cat("Model Equation Information \n")
  cat("\n")
  print(modeqns, quote = FALSE, right = FALSE, row.names = FALSE, print.gap=1)
  

}