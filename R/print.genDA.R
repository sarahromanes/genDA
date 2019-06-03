#' @title Print Method for genDA
#' @description Prints summarised information from the genDA function.
#'
#' @param object an object of class 'genDA' to print
#' @param ... Any other variables which will be ignored.
#' @export
#'
print.genDA <- function(object,...) {
  
  if (!inherits(object, "genDA"))  {
    stop("object not of class 'genDA'")
  }
  
 cat("                      genDA GLLVM fitting procedure                      \n")
 cat("\n")
 cat("\n")
 
 cat("Common Covariance Assumption: \n")
 if(names(object[1])=="call"){
   print(TRUE)
 } else {
   print(FALSE)
 }
 cat("\n")
 
 cat("Call: \n")
 if(names(object[1])=="call"){
  print(object$call)
 } else {
   print(object[[1]][[1]])
 }
 
 cat("\n")
 
 cat("log-likelihood: \n")
 if(names(object[1])=="call"){
   print(object$logL)
 } else {
   print(unlist(sapply(object, function(x) x[[2]])))
 }
 

 
}
