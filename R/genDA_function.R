genDA <- function(mY, mX, family, d=d){
  
  n <- nrow(mY)
  p <- ncol(mX)
  m <- ncol(mY)
  
  if(length(family)> 1 & m >1 & length(family)!=m){
    stop("Not all variables families accounted for, please check family vector")
  }
  if(length(family)==1 & m > 1){
    
    if(!(family %in% c("poisson","binomial", "gaussian", "log-normal")))
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
  }else{
    response_types <- family
    tmb_types <- rep(0, m)
    for(j in 1:m){
      
      if(!(response_types[j] %in% c("poisson","binomial", "gaussian", "log-normal")))
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
    }
  }
  
  vbeta0.hat <- c()
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
    vbeta0.hat[j] <- res$coef[1]
    mB.hat[,j] <- res$coef[-1]
  }
  
  # Set Ridge regression constants

  
  sigma2_beta0	<- 1E0
  vsigma2_beta    <- rep(1E1,m)
  vsigma2_lambda  <- rep(1.0E0,m)
  vsigma2_tau     <- rep(1.0E-2,n)
  
  vsigma2 <- c()
  vsigma2[1]           <- sigma2_beta0
  vsigma2[1+(1:m)]     <- vsigma2_beta
  vsigma2[1+m+(1:m)]   <- vsigma2_lambda
  vsigma2[1+m+m+(1:n)] <- vsigma2_tau
  
  
  vsigma_norm <- rep(0, m)
  for(j in 1:m){
    if(response_types[j]=="gaussian"|response_types[j]=="log-normal"){
      vsigma_norm[j] = sqrt(var(mY[,j]) *((n-1)/n))
    }
  }
  
  vtau.init <- rep(0,n)
  vbeta0.init <- vbeta0.hat
  mB.init <- mB.hat
  
  fit <-  genDA_start_values(mY, mX, family=family, d=d) 
  
  if (d > 0) {
    mL.init <- fit$mLambda
    if (d > 0)
      mL.init[upper.tri(mL.init)] <- 0
  }
  
  if(d > 0) lambda <- mL.init[lower.tri(mL.init,diag = TRUE)]
  
  mU.init <- fit$mU #extracts LV
  
  data <- list()
  data$mY <- mY
  data$mX <- mX
  data$vsigma2 <- vsigma2
  data$response_types <- tmb_types
  data$d <- d
  
  parameters <- list()
  parameters$log_vsigma_norm <- as.matrix(log(vsigma_norm))
  parameters$mB <- mB.init
  parameters$lambda <- lambda
  parameters$mU <- mU.init
  parameters$vtau <- as.matrix(vtau.init)
  parameters$vbeta0 <- as.matrix(vbeta0.init)
  
  
  obj = MakeADFun(data,parameters,DLL = "genDA_f", silent=TRUE, inner.control=list(mgcmax = 1e+200,maxit = 1000))
  opt1 = nlminb(obj$par,obj$fn, obj$gr,control=list(rel.tol=1.0E-6))
  
  param <- opt1$par
  li <- names(param)=="lambda"
  ui <- names(param)=="mU"
  bi <- names(param)=="mB"
  b0i <- names(param)=="vbeta0"
  ti <- names(param)=="vtau"
  vj <- names(param)=="log_vsigma_norm"
  
  val <- -1*opt1$objective
  
  if(d > 0){
    mU.hat<-(matrix(param[ui],n,d))
    theta <- matrix(0,m,d)
    if(m>1) {
      theta[lower.tri(theta,diag=TRUE)] <- param[li];
    } else {theta <- param[li]}
    
  }
  
  mB.hat <- t(matrix(nrow=m,ncol=p, param[bi]))
  colnames(mB.hat) <- colnames(mY)
  vtau.hat <- param[ti]
  vbeta0.hat <- opt1$par[b0i]
  
  inds <- which(response_types %in% c("gaussian", "log-normal"))
  if(length(inds)>0){
    vsigma_norm.hat <- rep(0,m)
    vsigma_norm.hat[inds] <- exp(opt1$par[vj][inds])
  }
  
  mEta.hat <- matrix(vtau.hat)%*%matrix(1,1,m) + matrix(1,n,1)%*%vbeta0.hat + mX%*%mB.hat + mU.hat%*%t(theta)
  
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
              vsigma_norm.hat = vsigma_norm,
              mU = mU.hat, 
              mEta = mEta.hat))
  
}