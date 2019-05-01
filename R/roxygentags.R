#' Roxygen commands
#'
#' This is a dummy function who's purpose is to hold the useDynLib roxygen tag.
#' This tag will populate the namespace with compiled c++ functions upon package install.
#'
#' @useDynLib genDA
#'
dummy <- function(){
  return(NULL)
}