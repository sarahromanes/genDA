#' Roxygen commands
#'
#' This is a dummy function who's purpose is to hold the useDynLib roxygen tag.
#' This tag will populate the namespace with compiled c++ functions upon package install.
#'
#' @useDynLib genDA_f
#' @useDynLib genDA_f_null_X
#' @useDynLib genDA_f_predict_approximation
#'
dummy <- function(){
  return(NULL)
}