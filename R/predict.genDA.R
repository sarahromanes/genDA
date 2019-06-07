#' @title Predict Method for genDA
#' @description Obtains new class predictions from a fitted genDA object.
#'
#' @export
#' @param object an object of class 'genDA'.
#' @param newdata A new data frame or matrix of response data. Must be numeric. Column response (family) types will be assumed to match that of the trained genDA object.
#' @param newX A new data frame or matrix of known covariate data.
#' @param ... not used.
#' 
#' @return An object including the following components:
#'
#'  \item{class }{a vector of predicted class labels}
#'  \item{prob_class }{a matrix of class membership probabilities}
#'  
#' @export


predict.genDA <- function(object, newdata, newX = NULL, ...){
  
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
  
  if(!is.null(newX)){
    newX <-  as.matrix(newX)
  }
  
  n_k = apply(.vec2mat(vc),2, sum)
  
  n_test <- nrow(newdata)
  m <- ncol(newdata)
  
  if(common.covariance){
    p <- object$side.list$p
    num.lv <- object$side.list$num.lv
    if(object$side.list$row.eff){vsigma2_tau <- rep(1.0E0,n_test)} else {vsigma2_tau <- rep(1.0E-4,n_test)}
    disp <-  object$side.list$disp
  } else {
    p <- object[[1]]$side.list$p
    num.lv <- object[[1]]$side.list$num.lv
    if(object[[1]]$side.list$row.eff){vsigma2_tau <- rep(1.0E0,n_test)} else {vsigma2_tau <- rep(1.0E-4,n_test)}
    disp <- object[[1]]$side.list$disp
  }

  row.init <- rep(0,n_test)
  
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
        
        if(object$side.list$row.eff){data$model_name = "genDA_f_predict"} else {data$model_name = "genDA_f_predict_null_row"}
        data$y <- t(as.matrix(newdata[i, ]))
        data$X <- t(as.matrix(c(class_test[k, ], newX[i, ])))
        if(object$side.list$row.eff){data$vsigma2_tau <- vsigma2_tau[i]}
        data$response_types <- object$side.list$tmb_types
        data$d <- num.lv
        data$mB <- object$params$Xcoef
        data$mL <- object$params$mL
        data$vbeta0 <- as.matrix(object$params$beta0)
        if(disp){data$vphi <- as.matrix(object$params$phi)} else{data$vphi <- as.matrix(rep(0, m))}
        
        parameters$mU <- t(as.matrix(mU.init[i,]))
        if(object$side.list$row.eff){parameters$vtau <- row.init[i]}
        
      } else {
        
        if(object[[k]]$side.list$row.eff){data$model_name = "genDA_f_predict"} else {data$model_name = "genDA_f_predict_null_row"}
        data$y <- t(as.matrix(newdata[i, ]))
        if(!is.null(newX)){data$X <- t(as.matrix(newX[i, ]))} else {data$X <- t(as.matrix(rep(0, 1)))}
        if(object[[k]]$side.list$row.eff){data$vsigma2_tau <- vsigma2_tau[i]}
        data$response_types <- object[[1]]$side.list$tmb_types
        data$d <- num.lv
        if(!is.null(newX)){data$mB <- object[[k]]$params$Xcoef} else {data$mB <- as.matrix(t(rep(0, m))) }
        data$mL <- object[[k]]$params$mL
        data$vbeta0 <- as.matrix(object[[k]]$params$beta0)
        if(disp){data$vphi <- as.matrix(object[[k]]$params$phi)} else{data$vphi <- as.matrix(rep(0, m))}
        
        parameters$mU <- t(as.matrix(mU.init[i,]))
        if(object[[k]]$side.list$row.eff){parameters$vtau <- row.init[i]}
        
      }
      
      obj <- MakeADFun(data,parameters,DLL = "genDA", silent=TRUE,  inner.control=list(mgcmax = 1e+200,maxit = 1000,tol10=0.01))
      opt <-  nlminb(obj$par,obj$fn, obj$gr,control=list(rel.tol=1.0E-6))
  
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
