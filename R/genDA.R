#' @title genDA: Multi-distributional Discriminant Analysis using Generalised Linear Latent Variable Modelling
#' @description Fits a genDA model for multivariate data through generalised linear latent variable modelling.
#'
#' @param y (n x m) matrix or data.frame of responses.
#' @param X matrix or data.frame of covariates
#' @param class numeric or factor of class information.
#' @param d  number of latent variables, d, in gllvm model. Non-negative integer, less than number of response variables (m). Defaults to 2.
#' @param family  distribution function for responses, to describe the distribution of each column. Columns can be of different family types. Family options are \code{poisson(link = "log")}, \code{"negative-binomial"} (with log link), \code{binomial(link = "logit")}, \code{"gaussian"}, and \code{"log-normal"}. Either a vector of family types matching the column length of the data can be provided, otherwise if a single family type is provided the algorithm will assume all columns match the single family type.
#' @param standard.errors logical. If \code{TRUE} (default is \code{FALSE}) standard errors for parameter estimates are calculated
#' 
#' @return An object of class "genDA" includes the following components:
#' 
#'  \item{call }{function call}
#'  \item{logL }{log likelihood}
#'  \item{lvs }{latent variables}
#'  \item{params}{list of parameters
#'  \itemize{
#'    \item{mL }{ coefficients related to latent variables}
#'    \item{vbeta0 }{ column specific intercepts}
#'    \item{mB }{ coefficients related to covariates X and/or class (X is combined with class if included)}
#'    \item{vtau }{ row-specific intercepts}
#'    \item{vphi }{ dispersion parameters \eqn{\phi} for negative binomial, Gaussian, or log-normal family columns. 0 for columns not belonging to these families}
#'    \item{mEta } calculated value of \eqn{\eta} with estimated coefficients above.
#'    }}
#'  \item{sd } standard errors of parameters
#'  \item{side.list }{list of additional parameters necessary to be passed onto to the predict.genDA method}

#'
#' @author Sarah Romanes <sarah.romanes@@sydney.edu.au>, John Ormerod
#' 
#' @export
#' @importFrom TMB MakeADFun
#' @importFrom mvabund manyglm
#' @importFrom stats glm nlminb var prcomp rnorm
#' @importFrom mvtnorm rmvnorm
#' @importFrom matrixStats colMeans2
#' 
genDA <- function(y, X = NULL, class = NULL, family, d=2, standard.errors = FALSE){
  
  y <- as.matrix(y)
  
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
    if(is.factor(class)){
      l <- levels(class)
      class <- as.character(class)
      class[which(class==l[1])] <-  0
      class[which(class==l[2])] <-  1
      class <- as.matrix(as.numeric(class))
    }
    
    if(!is.null(X)){
      X <- cbind(class, as.matrix(X))
    } else {
      X <- class
    }
    
    p <- ncol(X)
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
  
  vbeta0.hat <- c()
  
  if(!is.null(X)){
    mB.hat <- matrix(0,p,m)
    for (j in 1:m) {
      
      if (response_types[j]=="binomial") {
        res <- glm(y[,j]~X  ,family=binomial)
      }
      if (response_types[j]=="poisson") {
        res <- glm(y[,j]~X, family=poisson)
      }
      if(response_types[j]=="gaussian"){
        res <- lm(y[,j]~X)
      }
      if(response_types[j]=="log-normal"){
        res <- lm(log(y[,j])~X)
      }
      if(response_types[j]=="negative-binomial"){
        res <- mvabund::manyglm(y[,j]~X, family="negative.binomial", K=1)
      }
      vbeta0.hat[j] <- res$coef[1]
      mB.hat[,j] <- res$coef[-1]
    }
    
    # SET RIDGE REGRESSION CONSTANTS
    
    sigma2_beta0	<- 1E2
    vsigma2_beta    <- rep(1E1,m)
    vsigma2_lambda  <- rep(1.0E0,m)
    vsigma2_tau     <- rep(1.0E0,n)
    
    vsigma2 <- c()
    vsigma2[1]           <- sigma2_beta0
    vsigma2[1+(1:m)]     <- vsigma2_beta
    vsigma2[1+m+(1:m)]   <- vsigma2_lambda
    vsigma2[1+m+m+(1:n)] <- vsigma2_tau
    
  } else {
    for (j in 1:m) {
      
      if (response_types[j]=="binomial") {
        res <- glm(y[,j]~1  ,family=binomial)
      }
      if (response_types[j]=="poisson") {
        res <- glm(y[,j]~1, family=poisson)
      }
      if(response_types[j]=="gaussian"){
        res <- lm(y[,j]~1)
      }
      if(response_types[j]=="log-normal"){
        res <- lm(log(y[,j])~1)
      }
      if(response_types[j]=="negative-binomial"){
        res <- mvabund::manyglm(y[,j]~1, family="negative.binomial", K=1)
      }
      vbeta0.hat[j] <- res$coef[1]
  
    }
    
    # SET RIDGE REGRESSION CONSTANTS
    
    sigma2_beta0	<- 1E2
    vsigma2_lambda  <- rep(1.0E0,m)
    vsigma2_tau     <- rep(1.0E0,n)
    
    vsigma2 <- c()
    vsigma2[1]           <- sigma2_beta0
    vsigma2[1+(1:m)]   <- vsigma2_lambda
    vsigma2[1+m+(1:n)] <- vsigma2_tau
    
    }
  
  ## Initialise phi's ##
  
  vphi <- rep(0, m)
  label_phi <- vphi
  for(j in 1:m){
    if(response_types[j]=="gaussian"){
      vphi[j] <-  sqrt(var(y[,j]) *((n-1)/n))
      label_phi[j] <-  1 
    }
    if(response_types[j]=="log-normal"){
      vphi[j] <-  sqrt(var(log(y[,j])) *((n-1)/n))
      label_phi[j] <-  1
    }
    if(response_types[j]=="negative-binomial"){
      if(!is.null(X)){
        vphi[j] <-  mvabund::manyglm(y[,j]~X, family="negative.binomial", K=1)$phi
      } else {
        vphi[j] <-  mvabund::manyglm(y[,j]~1, family="negative.binomial", K=1)$phi
      }
      
      label_phi[j] <- 1
    }
  }

  
  fit <-  .genDA_start_values(y, X, family=family, d=d) 
  
  if (d > 0) {
    mL.init <- fit$mLambda
    if (d > 0)
      mL.init[upper.tri(mL.init)] <- 0
  }
  
  if(d > 0) lambda <- mL.init[lower.tri(mL.init,diag = TRUE)]
  
  mU.init <- fit$mU #extracts LV
  vtau.init <- rep(0,n)
  vbeta0.init <- vbeta0.hat
  if(!is.null(X)){mB.init <- mB.hat}
  
  data <- list()
  if(!is.null(X)){data$model_name="genDA_f"} else {data$model_name = "genDA_f_null_X"}
  data$y <- y
  if(!is.null(X)){data$X <- X}
  data$vsigma2 <- vsigma2
  data$response_types <- tmb_types
  data$d <- d
    
  parameters <- list()
  parameters$log_vphi <- as.matrix(log(vphi))
  if(!is.null(X)) {parameters$mB <- mB.init}
  parameters$lambda <- lambda
  parameters$mU <- mU.init
  parameters$vtau <- as.matrix(vtau.init)
  parameters$vbeta0 <- as.matrix(vbeta0.init)
    
  obj <-  MakeADFun(data,parameters,DLL = "genDA", silent=TRUE, inner.control=list(mgcmax = 1e+200,maxit = 1000))
  opt <-  nlminb(obj$par,obj$fn, obj$gr,control=list(rel.tol=1.0E-6))
  
  param <- opt$par
  lj <- names(param)=="lambda"
  ui <- names(param)=="mU"
  b0j <- names(param)=="vbeta0"
  ti <- names(param)=="vtau"
  vj <- names(param)=="log_vphi"
  if(!is.null(X)){bj <- names(param)=="mB"}
  
  val <- -1*opt$objective
  
  
  if(d > 0){
    mU.hat<-(matrix(param[ui],n,d))
    theta <- matrix(0,m,d)
    if(m>1) {
      theta[lower.tri(theta,diag=TRUE)] <- param[lj]
      rownames(theta) <- labels
    } else {theta <- param[lj]}
    colnames(theta)<-paste(rep("LV",ncol(theta)),c(1:ncol(theta)),sep="")
    colnames(mU.hat) <- colnames(theta)
  }
  
  if(!is.null(X)){
    mB.hat <- t(matrix(nrow=m,ncol=p, param[bj]))
    colnames(mB.hat) <- labels
  }
  
  vtau.hat <- param[ti]
  vbeta0.hat <- opt$par[b0j]; names(vbeta0.hat) <- labels
  
  
  vphi.hat <- rep(0,m)
  if(length(which(label_phi==1))>0){
    vphi.hat[which(label_phi==1)] <- exp(opt$par[vj][which(label_phi==1)])
  }
  names(vphi.hat) <- labels
  
  if(!is.null(X)){
    mEta.hat <- matrix(vtau.hat)%*%matrix(1,1,m) + matrix(1,n,1)%*%vbeta0.hat + X%*%mB.hat + mU.hat%*%t(theta)
  } else {
    mEta.hat <- matrix(vtau.hat)%*%matrix(1,1,m) + matrix(1,n,1)%*%vbeta0.hat  + mU.hat%*%t(theta)
  } 
  
  
  if(standard.errors==TRUE){
    
    sds <- obj$he()
    colnames(sds) <- names(param)
    rownames(sds) <- names(param)
    
    sd <- list()
    
    if(length(unique(response_types %in% c("negative-binomial", "gaussian", "log-normal")))==1){
      
      cov.mat <- solve(sds)
      ses <- sqrt(diag(abs(cov.mat)))
 
       if(all(response_types %in% c("negative-binomial", "gaussian", "log-normal"))){sd$vphi <- ses[names(ses)=="log_vphi"]*vphi.hat; names(sd$vphi) <- labels }
      
    } else {
  
      inds <- which(label_phi!=1)
      
      cov.mat <- solve(sds[-inds, -inds])
      ses <- sqrt(diag(abs(cov.mat)))
      
      sd$vphi <- ses[names(ses)=="log_vphi"]*vphi.hat[which(label_phi==1)]; names(sd$vphi) <- labels[-inds]
      
    }
    
    if(d > 0){
      sd$mU<-matrix(ses[names(ses)=="mU"],n,d)
      sd$mL <- matrix(0,m,d)
      if(m>1) {
        sd$mL[lower.tri(sd$mL,diag=TRUE)] <-  ses[names(ses)=="lambda"]
        rownames(sd$mL) <- labels
      } else {sd$mL <-  ses[names(ses)=="lambda"]}
      colnames(sd$mL)<-paste(rep("LV",ncol(theta)),c(1:ncol(theta)),sep="")
      colnames(sd$mU) <- colnames(sd$mL)
    }
    
    sd$vtau <- ses[names(ses)=="vtau"]; names(sd$vtau) <- labels.row
    sd$vbeta0 <- ses[names(ses)=="vbeta0"]; names(sd$vbeta0) <- labels 
    if(!is.null(X)){
      sd$mB <- t(matrix(nrow=m,ncol=p, ses[bj]))
      colnames(sd$mB) <- labels
    }
  }
  
  if(is.null(X)){
    mB.hat <-  NULL
    p <-  0
    vc <- NULL
  } else if (!is.null(X) & !is.null(class)){
    vc <- class
  } else {
    vc <- NULL
  }
  
  params <- list()
  params$mL <-t(theta)
  params$vbeta0 <- vbeta0.hat
  params$mB <- mB.hat
  params$vtau <- vtau.hat
  params$vphi <- vphi.hat
  params$mEta <- mEta.hat
  
  side.list <- list()
  side.list$p <- p
  side.list$tmb_types <- tmb_types
  side.list$vc <- vc
  side.list$n <- n
  side.list$d <- d
  
  object <- list()
  object$call <- match.call()
  object$logL <- val 
  object$lvs <- mU.hat
  object$params <- params
  if(standard.errors==TRUE) {object$sd <- sd} else {object$sd <-NULL} 
  object$side.list <- side.list
  
  class(object) <- "genDA"
  
  return(object)
  
}