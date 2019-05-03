#' @title Print Method for genDA
#' @description Prints summarised information from the genDA function.
#'
#' @param object an object of class 'genDA' to print
#' @param ... Any other variables which will be ignored.
#' @export
#'
print.multiDA <- function(object,...) {
  
  if (!inherits(object, "genDA"))  {
    stop("object not of class 'genDA'")
  }
  
 cat("Call: \n")
 print(object$call)
 
 cat("log-likelihood: \n")
 print(object$logL)
}
