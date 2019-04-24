genDA <- function(mY, mX, family, d=d){
  
  mY <- as.matrix(mY)
  
  if(!is.null(mX)){
    if(is.factor(mX)){
      l <- levels(mX)
      mX <- as.character(mX)
      mX[which(mX==l[1])] <-  0
      mX[which(mX==l[2])] <-  1
      mX <- as.matrix(as.numeric(mX))
    }
    
    p <- ncol(mX)
  }
  
  n <- nrow(mY)
  m <- ncol(mY)
  
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
        if(sum(mY[,j]<=0)>=1){
          stop(paste("Cannot model variable", j, "by Log-Normal - data must be strictly positive"))
        }
      }
      if(response_types[j]=="negative-binomial"){
        tmb_types[j] <- 5
      }
    }
  }
  
  vbeta0.hat <- c()
  
  if(!is.null(mX)){
    mB.hat <- matrix(0,p,m)
    for (j in 1:m) {
      
      if (response_types[j]=="binomial") {
        res <- glm(mY[,j]~mX  ,family=binomial)
      }
      if (response_types[j]=="poisson") {
        res <- glm(mY[,j]~mX, family=poisson)
      }
      if(response_types[j]=="gaussian"){
        res <- lm(mY[,j]~mX)
      }
      if(response_types[j]=="log-normal"){
        res <- lm(log(mY[,j])~mX)
      }
      if(response_types[j]=="negative-binomial"){
        #res <- MASS::glm.nb(mY[,j]~mX)
        res <- mvabund::manyglm(mY[,j]~mX, family="negative.binomial", K=1)
      }
      vbeta0.hat[j] <- res$coef[1]
      mB.hat[,j] <- res$coef[-1]
    }
    
    # SET RIDGE REGRESSION CONSTANTS
    
    sigma2_beta0	<- 1E2
    vsigma2_beta    <- rep(1E1,m)
    vsigma2_lambda  <- rep(1.0E0,m)
    vsigma2_tau     <- rep(1.0E-2,n)
    
    vsigma2 <- c()
    vsigma2[1]           <- sigma2_beta0
    vsigma2[1+(1:m)]     <- vsigma2_beta
    vsigma2[1+m+(1:m)]   <- vsigma2_lambda
    vsigma2[1+m+m+(1:n)] <- vsigma2_tau
    
  } else {
    for (j in 1:m) {
      
      if (response_types[j]=="binomial") {
        res <- glm(mY[,j]~1  ,family=binomial)
      }
      if (response_types[j]=="poisson") {
        res <- glm(mY[,j]~1, family=poisson)
      }
      if(response_types[j]=="gaussian"){
        res <- lm(mY[,j]~1)
      }
      if(response_types[j]=="log-normal"){
        res <- lm(log(mY[,j])~1)
      }
      if(response_types[j]=="negative-binomial"){
        res <- mvabund::manyglm(mY[,j]~1, family="negative.binomial", K=1)
      }
      vbeta0.hat[j] <- res$coef[1]
  
    }
    
    # SET RIDGE REGRESSION CONSTANTS
    
    sigma2_beta0	<- 1E2
    vsigma2_lambda  <- rep(1.0E0,m)
    vsigma2_tau     <- rep(1.0E-2,n)
    
    vsigma2 <- c()
    vsigma2[1]           <- sigma2_beta0
    vsigma2[1+(1:m)]   <- vsigma2_lambda
    vsigma2[1+m+(1:n)] <- vsigma2_tau
    
    }
  
  ## Initialise phi's ##
  
  vphi <- rep(0, m)
  for(j in 1:m){
    if(response_types[j]=="gaussian"){
      vphi[j] = sqrt(var(mY[,j]) *((n-1)/n))
    }
    if(response_types[j]=="log-normal"){
      vphi[j] = sqrt(var(log(mY[,j])) *((n-1)/n))
    }
    if(response_types[j]=="negative-binomial"){
      if(!is.null(mX)){
        vphi[j] = mvabund::manyglm(mY[,j]~mX, family="negative.binomial", K=1)$phi
      } else {
        vphi[j] = mvabund::manyglm(mY[,j]~1, family="negative.binomial", K=1)$phi
      }
    }
  }
  
  
  fit <-  genDA_start_values(mY, mX, family=family, d=d) 
  
  if (d > 0) {
    mL.init <- fit$mLambda
    if (d > 0)
      mL.init[upper.tri(mL.init)] <- 0
  }
  
  if(d > 0) lambda <- mL.init[lower.tri(mL.init,diag = TRUE)]
  
  mU.init <- fit$mU #extracts LV
  vtau.init <- rep(0,n)
  vbeta0.init <- vbeta0.hat
  if(!is.null(mX)){mB.init <- mB.hat}
  
  if(!is.null(mX)){
  
    data <- list()
    data$mY <- mY
    data$mX <- mX
    data$vsigma2 <- vsigma2
    data$response_types <- tmb_types
    data$d <- d
    
    parameters <- list()
    parameters$log_vphi <- as.matrix(log(vphi))
    parameters$mB <- mB.init
    parameters$lambda <- lambda
    parameters$mU <- mU.init
    parameters$vtau <- as.matrix(vtau.init)
    parameters$vbeta0 <- as.matrix(vbeta0.init)
    
    obj = MakeADFun(data,parameters,DLL = "genDA_f", silent=TRUE, inner.control=list(mgcmax = 1e+200,maxit = 1000))
  } else {
    
    data <- list()
    data$mY <- mY
    data$vsigma2 <- vsigma2
    data$response_types <- tmb_types
    data$d <- d
    
    parameters <- list()
    parameters$log_vphi <- as.matrix(log(vphi))
    parameters$lambda <- lambda
    parameters$mU <- mU.init
    parameters$vtau <- as.matrix(vtau.init)
    parameters$vbeta0 <- as.matrix(vbeta0.init)
    
    obj = MakeADFun(data,parameters,DLL = "genDA_f_null_X", silent=TRUE, inner.control=list(mgcmax = 1e+200,maxit = 1000))
    
  }
  
  opt = nlminb(obj$par,obj$fn, obj$gr,control=list(rel.tol=1.0E-6))
  
  param <- opt$par
  li <- names(param)=="lambda"
  ui <- names(param)=="mU"
  b0i <- names(param)=="vbeta0"
  ti <- names(param)=="vtau"
  vj <- names(param)=="log_vphi"
  if(!is.null(mX)){bi <- names(param)=="mB"}
  
  val <- -1*opt$objective
  
  if(d > 0){
    mU.hat<-(matrix(param[ui],n,d))
    theta <- matrix(0,m,d)
    if(m>1) {
      theta[lower.tri(theta,diag=TRUE)] <- param[li];
    } else {theta <- param[li]}
    
  }
  
  if(!is.null(mX)){
    mB.hat <- t(matrix(nrow=m,ncol=p, param[bi]))
    colnames(mB.hat) <- colnames(mY)
  }
  
  vtau.hat <- param[ti]
  vbeta0.hat <- opt$par[b0i]
  
  inds <- which(response_types %in% c("gaussian", "log-normal", "negative-binomial"))
  if(length(inds)>0){
    vphi.hat <- rep(0,m)
    vphi.hat[inds] <- exp(opt$par[vj][inds])
  }
  
  if(!is.null(mX)){
    mEta.hat <- matrix(vtau.hat)%*%matrix(1,1,m) + matrix(1,n,1)%*%vbeta0.hat + mX%*%mB.hat + mU.hat%*%t(theta)
  } else {
    mEta.hat <- matrix(vtau.hat)%*%matrix(1,1,m) + matrix(1,n,1)%*%vbeta0.hat  + mU.hat%*%t(theta)
  } 
  
  if(is.null(mX)){
    mB.hat <-  NULL
    p <-  0
    vc <- NULL
  }
  
  return(list(mB = mB.hat, 
              mL = t(theta), 
              p = p, 
              family = family, 
              tmb_types = tmb_types, 
              d = d,
              ll = val,
              vc = mX[,1],
              n = n,
              vtau = vtau.hat,
              vbeta0 = vbeta0.hat,
              vphi = vphi.hat,
              mU = mU.hat, 
              mEta = mEta.hat))
  
}