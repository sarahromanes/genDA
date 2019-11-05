
####################### LDA FIT (common covariance structure) ###################################

.genDA_fit_LDA <- function(y, X, class, num.lv, tmb_types, response_types, standard.errors, call, labels, labels.row, family){
  
  n <- nrow(y)
  m <- ncol(y)
  
  if(!is.null(X)){
    p <- ncol(X)
    mB.hat <- matrix(0,p,m)
    for (j in 1:m) {
      
      if (response_types[j]=="binomial") {
        res <- glm(y[,j]~X-1  ,family=binomial)
      }
      if (response_types[j]=="poisson") {
        res <- glm(y[,j]~X-1, family=poisson)
      }
      if(response_types[j]=="gaussian"){
        res <- lm(y[,j]~X-1)
      }
      if(response_types[j]=="log-normal"){
        res <- lm(log(y[,j])~X-1)
      }
      if(response_types[j]=="negative-binomial"){
        res <- mvabund::manyglm(y[,j]~X-1, family="negative.binomial", K=1)
      }
      mB.hat[,j] <- res$coef
    }
    
    # SET RIDGE REGRESSION CONSTANTS
    
    vsigma2_beta    <- rep(1E1,m)
    vsigma2_lambda  <- rep(1.0E1,m)
    
  } else {
    
    vbeta0.hat <- c()
    
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
    
    sigma2_beta0	<- 1E0
    vsigma2_lambda  <- rep(1.0E0,m)

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
  
  disp <- !all(label_phi==0)
  if(disp){
    vphi <- vphi[label_phi==1]
    vphi_inds <- rep(0,m)
    vphi_inds[which(label_phi==1)]= 1:length(which(label_phi==1))
    vphi_inds <- vphi_inds - 1
  }
  
  
  fit <-  .genDA_start_values(y, X, family=family, num.lv=num.lv) 
  
  if (num.lv > 0) {
    mL.init <- fit$mLambda
    if (num.lv > 0)
      mL.init[upper.tri(mL.init)] <- 0
  }
  
  if(num.lv > 0) lambda <- mL.init[lower.tri(mL.init,diag = TRUE)]
  
  mU.init <- fit$mU #extracts LV
  if(is.null(X)){vbeta0.init <- vbeta0.hat}
  if(!is.null(X)){mB.init <- mB.hat}
  
  # ----------------- FIT TMB ----------------------- # 
  
  if(!is.null(X) & disp){
    model_name = "genDA_f"
  } else if(!is.null(X) & !disp){
    model_name = "genDA_f_null_dis"
  } else if(is.null(X) & disp){
    model_name = "genDA_f_null_X"
  } else if(is.null(X) & !disp){
    model_name = "genDA_f_null_X_dis"
  } 
  
  data <- list()
  data$model_name <- model_name
  data$y <- y
  if(!is.null(X)){data$X <- X}
  if(is.null(X)){data$sigma2_beta0 <- sigma2_beta0}
  if(!is.null(X)){data$vsigma2_beta <- vsigma2_beta}
  data$vsigma2_lambda <- vsigma2_lambda
  data$response_types <- tmb_types
  data$d <- num.lv
  if(!is.null(X)){data$p <- nrow(mB.init)}
  if(disp){data$vphi_inds <- vphi_inds}
  
  parameters <- list()
  if(disp){parameters$log_vphi <- as.matrix(log(vphi))}
  #if(!is.null(X)) {parameters$mB <- mB.init}
  if(!is.null(X)){parameters$beta <- as.numeric(t(mB.init))}
  parameters$lambda <- lambda
  parameters$mU <- mU.init
  if(is.null(X)){parameters$vbeta0 <- as.matrix(vbeta0.init)}
  
  
  obj <-  MakeADFun(data,parameters,DLL = "genDA", silent=TRUE, inner.control=list(mgcmax = 1e+200,maxit = 1000))
  opt <-  nlminb(obj$par,obj$fn, obj$gr,control=list(rel.tol=1.0E-6))
  
  param <- opt$par
  lj <- names(param)=="lambda"
  ui <- names(param)=="mU"
  if(is.null(X)){b0j <- names(param)=="vbeta0"}
  if(disp){vj <- names(param)=="log_vphi"}
  if(!is.null(X)){bj <- names(param)=="beta"}
  
  val <- -1*opt$objective
  
  
  if(num.lv > 0){
    mU.hat<-(matrix(param[ui],n,num.lv))
    theta <- matrix(0,m,num.lv)
    if(m>1) {
      theta[lower.tri(theta,diag=TRUE)] <- param[lj]
      rownames(theta) <- labels
    } else {theta <- param[lj]}
    colnames(theta)<-paste(rep("LV",ncol(theta)),c(1:ncol(theta)),sep="")
    colnames(mU.hat) <- colnames(theta)
  }
  
  mB.hat <-  NULL
  if(!is.null(X)){
    beta.hat <- param[bj]
    mB.hat <- matrix(beta.hat, p, m, byrow =T)
    colnames(mB.hat) <- labels
    rownames(mB.hat) <- colnames(X)
    vbeta0.hat <- mB.hat[1,] ; names(vbeta0.hat) <- labels
    mB.hat <- as.matrix(mB.hat[-1,])
    p <- ncol(mB.hat)
  }
  
  if(is.null(X)){vbeta0.hat <- opt$par[b0j]; names(vbeta0.hat) <- labels}
  
  vphi.hat <- NULL
  if(disp){
    vphi.hat <- rep(0,m)
    for(j in 1:m){
      if(label_phi[j]!=1){next}
 
        vphi.hat[which(label_phi==1)] <- exp(opt$par[vj])
      
    }
    names(vphi.hat) <- labels
  }

  if(standard.errors){
    
    sd <- list()
    
    sds <- optimHess(obj$par,obj$fn, obj$gr)
    cov.mat <- ginv(sds)
    ses <- sqrt(diag(abs(cov.mat)))
    names(ses) <- names(param)
    
    
    if(disp){
      for(j in 1:m){
        if(label_phi[j]!=1){
          next
        }
        
          sd$phi <- ses[names(ses)=="log_vphi"]*vphi.hat/(1 + vphi.hat); names(sd$phi) <- labels
        
      }
    } 
    
    if(num.lv > 0){
      sd$mU<-matrix(ses[names(ses)=="mU"],n,num.lv)
      sd$mL <- matrix(0,m,num.lv)
      if(m>1) {
        sd$mL[lower.tri(sd$mL,diag=TRUE)] <-  ses[names(ses)=="lambda"]
        rownames(sd$mL) <- labels
      } else {sd$mL <-  ses[names(ses)=="lambda"]}
      colnames(sd$mL)<-paste(rep("LV",ncol(theta)),c(1:ncol(theta)),sep="")
      colnames(sd$mU) <- colnames(sd$mL)
    }
    
    if(is.null(X)){sd$beta0 <- ses[names(ses)=="vbeta0"]; names(sd$beta0) <- labels}
    
    if(!is.null(X)){
      sd.beta <- matrix(ses[names(ses)=="beta"], p+1, m, byrow = T)
      sd$vbeta0 <- sd.beta[1, ]; names(sd$vbeta0) <- labels
      sd$Xcoef <- t(as.matrix(sd.beta[-1,]))
      colnames(sd$Xcoef) <- labels
    }
    
  }
  
  if(is.null(X)){p <-  0}
  if (!is.null(class)){
    vc <- class
  } else if (is.null(class)) {
    vc <- NULL
  }
  
  params <- list()
  params$mL <-t(theta)
  params$beta0 <- vbeta0.hat
  params$Xcoef <- mB.hat
  params$phi <- vphi.hat
  
  side.list <- list()
  side.list$p <- p
  side.list$tmb_types <- tmb_types
  side.list$vc <- vc
  side.list$n <- n
  side.list$num.lv <- num.lv
  side.list$disp <- disp
  
  object <- list()
  object$call <- call
  object$logL <- val 
  object$lvs <- mU.hat
  object$params <- params
  if(standard.errors==TRUE) {object$sd <- sd} else {object$sd <-NULL} 
  object$side.list <- side.list
  
  class(object) <- "genDA"
  
  
  return(object)
}


####################### QDA FIT (seperate covariance structure) #################################

.genDA_fit_QDA <- function(y, class, num.lv, tmb_types, response_types, standard.errors, call, labels, labels.row, family){
  
  K <- length(levels(class))
  
  object <- list()
 
  for(k in 1:K){
    y_sub <- y[which(class==levels(class)[k]), ]
    
    object[[k]] <- quiet(.genDA_fit_LDA(y =y_sub, X =NULL, class = class, num.lv = num.lv,tmb_types = tmb_types, response_types =  response_types, standard.errors = standard.errors, call = call, labels =labels, labels.row =labels.row[which(class==levels(class)[k])], family=family))
    
  }
  
  names(object) <- levels(class)
  
  class(object) <- "genDA"
  
  return(object)
}

