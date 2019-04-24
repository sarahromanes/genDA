genDA_start_values<- function(mY, mX = NULL, family, d = 0) {
  
  set.seed(1)

  N <-n <- nrow(mY); 
  m <- ncol(mY); 
  mY <- as.matrix(mY)
  
  num.X <- 0; 
  if(!is.null(mX)) num.X <- dim(mX)[2]
  
  if(!is.numeric(mY))
    stop("mY must be numeric")
  
  if(d > 0) {
    mU <- mvtnorm::rmvnorm(N, rep(0, d));
  }
  
  if(d == 0) { mU <- NULL }

  if(length(family)==1 & m > 1){
    
    if(!(family %in% c("poisson","binomial", "gaussian", "log-normal", "negative-binomial")))
      stop("Input family not supported")
    
    if(family=="log-normal"){
      family="gaussian"
      mY <- log(mY)
    }
    
    if(family=="negative-binomial"){
      family="negative.binomial"
    }
    
    if(!is.null(mX)){
      if(family=="gaussian"){
        fit.mva <- mvabund::manylm(mY~mX)
      } else {
        fit.mva <- mvabund::manyglm(mY ~ mX, family = family, K = 1)
      }
    }
    if(is.null(mX)) {
      if(family=="gaussian"){
        fit.mva <- mvabund::manylm(mY ~ 1)
      } else {
        fit.mva <- mvabund::manyglm(mY ~ 1, family = family, K = 1) 
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
        mY[,j] <- log(mY[,j])
      }
      if(response_types[j]=="negative-binomial"){
        response_types[j]="negative.binomial"
      }
      
      if(!is.null(mX)){
        if(response_types[j]=="gaussian"){
          fit.mva <- mvabund::manylm(mY[,j]~mX)
        } else {
          fit.mva <- mvabund::manyglm(mY[,j] ~ mX, family = response_types[j], K = 1)
        }
      }
      if(is.null(mX)) {
        if(response_types[j]=="gaussian"){
          fit.mva <- mvabund::manylm(mY[,j] ~ 1)
        } else {
          fit.mva <- mvabund::manyglm(mY[,j] ~ 1, family = response_types[j], K = 1) 
        }
      }
      resi_j <- residuals(fit.mva); resi_j[is.infinite(resi_j)] <- 0; resi_j[is.nan(resi_j)] <- 0
      resi <- cbind(resi, resi_j)
    }
   
  }
  
  mLambda <- NULL
  if(d > 0){
    if(m > 2 && n > 2){
      
      fa <-  EM.FA(resi, d)
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
