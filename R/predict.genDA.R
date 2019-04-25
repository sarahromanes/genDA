#' @title Predict Method for genDA
#' @description Obtains new class predictions from a fitted genDA object.
#'
#' @param object an object of class 'genDA'.
#' @param new.y A new data frame or matrix of response data. Must be numeric. Column response (family) types will be assumed to match that of the trained genDA object.
#' @param newX A new data frame or matrix of known covariate data.
#' @param prior_beta A vector of parameters for the beta prior for class prediction. Defaults to c(1,1) (Uninformative prior)
#' @param ... not used.
#' 
#' @return An object including the following components:
#'
#'  \item{class }{a vector of predicted class labels}
#'  \item{prob_class }{a matrix of class membership probabilities}
#'  
#' @author Sarah Romanes <sarah.romanes@@sydney.edu.au>
#' 
#' @aliases predict predict.genDA
#' @method predict genDA
#' @export
#' @export predict.genDA

predict.genDA <- function(object, new.y, newX = NULL, prior_beta = c(1,1), ...){
  
  new.y <- as.matrix(new.y)
  newX <-  as.matrix(newX)
  
  if (!inherits(object, "genDA"))  {
    stop("object not of class 'genDA'")
  }
  
  vc <- object$predict.params$vc
  if(is.null(vc)){
    stop("GLLVM not trained with class responses included. Try fitting genDA again with class information captured 'class'.")
  }
  
  n_test <- nrow(new.y)
  p <- object$predict.params$p
  m <- ncol(new.y)
  d <- object$predict.params$d
  n_train <- object$predict.params$n
  
  A <- prior_beta[1]
  B <- prior_beta[2]

  vsigma2_tau     <- rep(1.0E-4,n_test)
  vtau.init <- rep(0,n_test)
  
  mU.init <- matrix(rnorm(n_test*d), nrow=n_test, ncol=d) # get same results if I use start.values.gllvm.TMB or random. 
  
  class <- rep(0, n_test)
  prob_class <- matrix(0,nrow = n_test, ncol=2)
  colnames(prob_class) <- c("Class 0", "Class 1")
  
  for(i in 1:n_test){
    
    data_1 <- list()
    data_1$y <- t(as.matrix(new.y[i, ]))
    data_1$X <- t(as.matrix(c(1, newX[i, ]))) #setting vc_test[i] = 1
    data_1$vsigma2_tau <- vsigma2_tau[i]
    data_1$response_types <- object$predict.values$tmb_types
    data_1$d <- d
    data_1$mB <- object$params$mB
    data_1$mL <- object$params$mL
    data_1$vbeta0 <- as.matrix(object$params$vbeta0)
    data_1$vphi <- as.matrix(object$params$vphi)
    
    data_0 <- data_1
    data_0$X <-  t(as.matrix(c(0, newX[i, ])))
        
    parameters <- list()
    parameters$mU <- t(as.matrix(mU.init[i,]))
    parameters$vtau <- vtau.init[i]
    
    
    obj_1 <- MakeADFun(data_1,parameters,DLL = "genDA_f_predict_approximation", silent=TRUE)
    obj_0 <-  MakeADFun(data_0,parameters,DLL = "genDA_f_predict_approximation", silent=TRUE)
    opt_1 <-  nlminb(obj_1$par,obj_1$fn, obj_1$gr,control=list(rel.tol=1.0E-6))
    opt_0 <- nlminb(obj_0$par,obj_0$fn, obj_0$gr,control=list(rel.tol=1.0E-6))
    
    ll_1 <- -1*(opt_1$objective) +  lbeta(A + sum(object$predict.values$vc) + 1, B + n_train - sum(object$predict.values$vc) -1) 
    ll_0 <- -1*(opt_0$objective) +  lbeta(A + sum(object$predict.values$vc), B + n_train - sum(object$predict.values$vc)) 
    
    if(ll_1 > ll_0 ){
      class[i] = 1
    }
    
    prob_class[i, 1] <- exp(ll_0)/(exp(ll_0)+exp(ll_1))
    prob_class[i, 2] <- 1- prob_class[i,1]
    
  }
  
  return(list(class = class, 
              prob_class = round(prob_class,4)))
  
}

#' @export predict
predict <- function(object, ...)
{
  UseMethod(generic = "predict")
} 