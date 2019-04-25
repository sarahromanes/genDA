## Function to perform Factor Analysis by EM method ##

.EM.FA <- function(y, k){
  
  n <- nrow(y)
  p <- ncol(y)
  
  pc <- prcomp(y)
  mD <- diag(p)
  mZ <- pc$rotation[,1:k]
  
  
  # Our estimate of vmu is calculated once out the iteration loop.
  
  vmu <- matrixStats::colMeans2(y)
  mMu <- matrix(rep(vmu, n), nrow=n, ncol=p, byrow=T)
  
  min.err <- 1.0E-8
  
  crit <- 1
  
  
  while(crit > min.err){
    
    mS <- solve(crossprod(mZ,solve(mD))%*%mZ + diag(k))
    mM <- tcrossprod(mS,as.matrix(mZ))%*%tcrossprod(solve(mD),y-mMu) # note, we have a matrix of m, each column for each observation
    
    eU <- mM
    
    eUU.sum <- matrix(0, nrow=k, ncol=k)
    for(j in 1:n){
      eUU.sum <-eUU.sum +  mS +tcrossprod(eU[,j],eU[,j])
    }
    
    eU.sum <- matrix(0,nrow=k, ncol=p)
    for(j in 1:n){
      eU.sum <-eU.sum+ tcrossprod(eU[,j], y[j,]-vmu)  
      
    }
    
    mZ.new <- crossprod(eU.sum,solve(eUU.sum))
    mD.new <- diag(diag(cov(y)*((n-1)/n) - (mZ.new%*%eU.sum/n))) #note I have corrected the sample cov matrix to data cov matrix
    
    crit <- mean((mD - mD.new)^2)
    
    mD <- mD.new
    mZ <- mZ.new
    
    
  }
  
  
  return(list(mD=mD, mZ=mZ, mM=mM))
}


## Function to initialise genDA starting values ##

.genDA_start_values<- function(y, X = NULL, family, d = 2) {
  
  set.seed(1)

  N <-n <- nrow(y); 
  m <- ncol(y); 
  y <- as.matrix(y)
  
  num.X <- 0; 
  if(!is.null(X)) num.X <- dim(X)[2]
  
  if(!is.numeric(y))
    stop("y must be numeric")
  
  if(d > 0) {
    mU <- mvtnorm::rmvnorm(N, rep(0, d));
  }
  
  if(d == 0) { mU <- NULL }

  if(length(family)==1 & m > 1){
    
    if(!(family %in% c("poisson","binomial", "gaussian", "log-normal", "negative-binomial")))
      stop("Input family not supported")
    
    if(family=="log-normal"){
      family="gaussian"
      y <- log(y)
    }
    
    if(family=="negative-binomial"){
      family="negative.binomial"
    }
    
    if(!is.null(X)){
      if(family=="gaussian"){
        fit.mva <- mvabund::manylm(y~X)
      } else {
        fit.mva <- mvabund::manyglm(y ~ X, family = family, K = 1)
      }
    }
    if(is.null(X)) {
      if(family=="gaussian"){
        fit.mva <- mvabund::manylm(y ~ 1)
      } else {
        fit.mva <- mvabund::manyglm(y ~ 1, family = family, K = 1) 
      }
    }
    
    resi <- residuals(fit.mva); resi[is.infinite(resi)] <- 0; resi[is.nan(resi)] <- 0
    
  }else{
    response_types <- family
    resi <- c()
    coef <- c()
    for(j in 1:m){
      if(response_types[j]=="log-normal"){
        response_types[j]="gaussian"
        y[,j] <- log(y[,j])
      }
      if(response_types[j]=="negative-binomial"){
        response_types[j]="negative.binomial"
      }
      
      if(!is.null(X)){
        if(response_types[j]=="gaussian"){
          fit.mva <- mvabund::manylm(y[,j]~X)
        } else {
          fit.mva <- mvabund::manyglm(y[,j] ~ X, family = response_types[j], K = 1)
        }
      }
      if(is.null(X)) {
        if(response_types[j]=="gaussian"){
          fit.mva <- mvabund::manylm(y[,j] ~ 1)
        } else {
          fit.mva <- mvabund::manyglm(y[,j] ~ 1, family = response_types[j], K = 1) 
        }
      }
      resi_j <- residuals(fit.mva); resi_j[is.infinite(resi_j)] <- 0; resi_j[is.nan(resi_j)] <- 0
      resi <- cbind(resi, resi_j)
    }
   
  }
  
  mLambda <- NULL
  if(d > 0){
    if(m > 2 && n > 2){
      
      fa <-  .EM.FA(resi, d)
      mLambda <- (fa$mZ)
      mU <- t(fa$mM)
      
    } else {
      mLambda <- matrix(1,m,d)
      mLambda[upper.tri(mLambda)]=0
      mU <- matrix(0,n,d)
    }
  }
  params <- mLambda
  
  if(d > 0) {
    mU <- mU+mvtnorm::rmvnorm(n, rep(0, d),diag(d))
  }

  out <- list(mLambda=mLambda)
  
  if(d > 0) out$mU <- mU
  
  return(out)
}


