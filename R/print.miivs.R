#' @method print miivs 
#' @export
print.miivs <- function(x,...){
  modeqns   <- x$df
  cat("Model Equation Information \n")
  cat("\n")
  print(modeqns, quote = FALSE, right = FALSE, row.names = FALSE, print.gap=1)
}