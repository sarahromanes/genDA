#' @title genDA: Multi-distributional Discriminant Analysis using Generalised Linear Latent Variable Modelling
#' @description Fits a genDA model for multivariate data through generalised linear latent variable modelling.
#'
#' @param y (n x m) matrix or data.frame of responses.
#' @param class a factor of class information.
#' @param num.lv  number of latent variables, d, in gllvm model. Non-negative integer, less than number of response variables (m). Defaults to 2.
#' @param family  distribution function for responses, to describe the distribution of each column. Columns can be of different family types. Family options are \code{"poisson"} (with log link), \code{"negative-binomial"} (with log link), \code{"binomial"} (with \code{logit} link), \code{"gaussian"}, and \code{"log-normal"}. Either a vector of family types matching the column length of the data can be provided, otherwise if a single family type is provided the algorithm will assume all columns match the single family type.
#' @param standard.errors logical. If \code{TRUE} (default is \code{FALSE}) standard errors for parameter estimates are calculated.
#' @param common.covariance logical. Default \code{TRUE}. Specifies whether different covariance structures should be fit for each class, if \code{class} is non-null. If class is not provided, then only one covariance structure can be fit. See vignette for more details.
#' 
#' @return An object of class "genDA" includes the following components below. If \code{common.covariance} is set to \code{FALSE}, this information will be contained in seperate lists for each of the models fit for the seperate classes.
#' 
#'  \item{call }{function call}
#'  \item{logL }{log likelihood}
#'  \item{lvs }{latent variables}
#'  \item{params}{list of parameters
#'  \itemize{
#'    \item{mL }{ coefficients related to latent variables}
#'    \item{beta0 }{ column specific intercepts}
#'    \item{Xcoef }{ coefficients related to covariates X and/or class (X is combined with class if included)}
#'    \item{row }{ row-specific intercepts}
#'    \item{phi }{ dispersion parameters \eqn{\phi} for negative binomial, Gaussian, or log-normal family columns. 0 for columns not belonging to these families}
#'    }}
#'  \item{sd } standard errors of parameters
#'  \item{side.list }{list of additional parameters necessary to be passed onto to the predict.genDA method}
#'
#' @author Sarah Romanes <sarah.romanes@@sydney.edu.au>, John Ormerod
#' 
#' @export
#' @importFrom TMB MakeADFun
#' @importFrom mvabund manyglm
#' @importFrom stats glm nlminb var prcomp rnorm optimHess
#' @importFrom mvtnorm rmvnorm
#' @importFrom matrixStats colMeans2
#' @importFrom MASS ginv
#' @importFrom purrr map_chr
 
genDA <- function(y, class = NULL, family, num.lv=2, standard.errors = FALSE, common.covariance = TRUE){
  
  y <- as.matrix(y)
  
  if(!common.covariance & is.null(class)){
    warning("Cannot fit seperate covariances if no classes are specified. Will run with common covariance")
    common.covariance = TRUE
  }
  
  
  if(is.null(colnames(y))){
    colnames(y)<-paste("V",1:ncol(y), sep="")
  }else{
    colnames(y)<-make.unique(colnames(y))
  }
  
  if(is.null(rownames(y))){
    rownames(y)<-paste("R",1:nrow(y), sep="")
  }else{
    rownames(y)<-make.unique(rownames(y))
  }
  
  labels <- colnames(y)
  labels.row <- rownames(y)
  
  if(!is.null(class)){
    
    K = length(unique(class))
    classF <- as.factor(class)
    class <- .vec2mat(classF)
    colnames(class) <- levels(classF)
    class[,1] <- 1 # ensure columns are linearly independent
    
    if(common.covariance){
      X <- as.matrix(class)
      p <- ncol(X)
    }
      
  } else {
    classF <- NULL
    X <- NULL
  }

  
  n <- nrow(y)
  m <- ncol(y)
  
  if(length(family)> 1 & m >1 & length(family)!=m){
    stop("Not all variables families accounted for, please check family vector")
  }
  if(length(family)==1 & m > 1){
    
    if(!(family %in% c("poisson","binomial", "gaussian", "log-normal", "negative-binomial")))
      stop("Input family not supported")
    
    response_types <- rep(family, each=m)
    
    if(family=="binomial"){
      tmb_types <- rep(1,m)
    }
    if(family=="poisson"){
      tmb_types <- rep(2,m)
    }
    if(family=="gaussian"){
      tmb_types <- rep(3,m)
    }
    if(family=="log-normal"){
      tmb_types <- rep(4,m)
    }
    if(family=="negative-binomial"){
      tmb_types <- rep(5,m)
    }
  }else{
    response_types <- family
    tmb_types <- rep(0, m)
    for(j in 1:m){
      
      if(!(response_types[j] %in% c("poisson","binomial", "gaussian", "log-normal", "negative-binomial")))
        stop(paste("Variable", j, "- input family :",response_types[j],"not supported"))
      
      if(response_types[j]=="binomial"){
        tmb_types[j] <- 1
      }
      if(response_types[j]=="poisson"){
        tmb_types[j] <- 2
      }
      if(response_types[j]=="gaussian"){
        tmb_types[j] <- 3
      }
      if(response_types[j]=="log-normal"){
        tmb_types[j] <- 4
        if(sum(y[,j]<=0)>=1){
          stop(paste("Cannot model variable", j, "by Log-Normal - data must be strictly positive"))
        }
      }
      if(response_types[j]=="negative-binomial"){
        tmb_types[j] <- 5
      }
    }
  }
  
 call <- match.call()
  
  if(common.covariance){
    return(suppressWarnings(.genDA_fit_LDA(y = y, X = X, class = classF, num.lv = num.lv, tmb_types = tmb_types, response_types =response_types, standard.errors = standard.errors, call = call, labels = labels, labels.row = labels.row, family = family)))
  } else {
    return(suppressWarnings(.genDA_fit_QDA(y = y, class = classF, num.lv = num.lv, tmb_types = tmb_types, response_types=response_types, standard.errors= standard.errors, call = call, labels = labels, labels.row = labels.row, family = family)))
  }

}