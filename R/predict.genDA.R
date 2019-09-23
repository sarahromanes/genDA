#' @title Predict Method for genDA
#' @description Obtains new class predictions from a fitted genDA object.
#'
#' @export
#' @param object an object of class 'genDA'.
#' @param newdata A new data frame or matrix of response data. Must be numeric. Column response (family) types will be assumed to match that of the trained genDA object.
#' @param ... not used.
#' 
#' @return An object including the following components:
#'
#'  \item{class }{a vector of predicted class labels}
#'  \item{prob_class }{a matrix of class membership probabilities}
#'  
#' @export


predict.genDA <- function(object, newdata, ...){
  
  if (!inherits(object, "genDA"))  {
    stop("object not of class 'genDA'")
  }
  
  if(names(object[1])=="call"){
    common.covariance = TRUE
  } else {
    common.covariance = FALSE
    class_names = names(object)
    vc <- object[[1]]$side.list$vc
  }
  
  if(common.covariance){
    vc <- object$side.list$vc
    if(is.null(vc)){
      stop("GLLVM not trained with class responses included. Try fitting genDA again with class information captured in 'class'.")
    }
    class_names <- levels(object$side.list$vc)
  }
  
  newdata <- as.matrix(newdata)
  
  n_k = apply(.vec2mat(vc),2, sum)
  
  n_test <- nrow(newdata)
  m <- ncol(newdata)
  
  if(common.covariance){
    p <- object$side.list$p
    num.lv <- object$side.list$num.lv
    disp <-  object$side.list$disp
  } else {
    p <- object[[1]]$side.list$p
    num.lv <- object[[1]]$side.list$num.lv
    disp <- object[[1]]$side.list$disp
  }

  mU.init <- matrix(rnorm(n_test*num.lv), nrow=n_test, ncol=num.lv) 
  
  K <- length(class_names)
  ll <- matrix(0,nrow = n_test, ncol=K)
  colnames(ll) <- class_names
  
  class_test <-  as.matrix(diag(K)[, -1])
  
  
  for(i in 1:n_test){
  
    
    for(k in 1:K){
      
      data <- list()
      parameters <- list()
      
      if(common.covariance){
        
        data$model_name = "genDA_f_predict"
        data$y <- t(as.matrix(newdata[i, ]))
        data$X <- t(as.matrix(class_test[k, ]))
        data$response_types <- object$side.list$tmb_types
        data$d <- num.lv
        data$mB <- as.matrix(object$params$Xcoef)
        data$mL <- object$params$mL
        data$vbeta0 <- as.matrix(object$params$beta0)
        if(disp){data$vphi <- as.matrix(object$params$phi)} else{data$vphi <- as.matrix(rep(0, m))}
        
        parameters$mU <- t(as.matrix(mU.init[i,]))
        
      } else {
        
        data$model_name = "genDA_f_predict"
        data$y <- t(as.matrix(newdata[i, ]))
        data$X <- t(as.matrix(rep(0, 1)))
        data$response_types <- object[[1]]$side.list$tmb_types
        data$d <- num.lv
        data$mB <- as.matrix(t(rep(0, m))) 
        data$mL <- object[[k]]$params$mL
        data$vbeta0 <- as.matrix(object[[k]]$params$beta0)
        if(disp){data$vphi <- as.matrix(object[[k]]$params$phi)} else{data$vphi <- as.matrix(rep(0, m))}
        
        parameters$mU <- t(as.matrix(mU.init[i,]))

        
      }
      
      obj <- suppressWarnings(MakeADFun(data,parameters,DLL = "genDA", silent=TRUE,  inner.control=list(mgcmax = 1e+200,maxit = 1000,tol10=0.01)))
      opt <-  suppressWarnings(nlminb(obj$par,obj$fn, obj$gr,control=list(rel.tol=1.0E-6)))
  
      n_k[k] <- n_k[k] + 1
      
      ll[i,k]<- -1*(opt$objective) +  sum(lgamma(n_k)) - lgamma(sum(n_k))
      
      n_k[k] <- n_k[k] -1
    }
    
  }
  
  prob_class <- .logMatToGamma(ll)
  class <-  apply(prob_class, 1, which.max)
  class <-  as.factor(purrr::map_chr(class, .num.2.fac,vc))
  
  return(list(class = class, 
              prob_class = round(prob_class,4)))
  
}
