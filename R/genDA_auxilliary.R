
########################################################################################

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
    
    mS <- tryCatch(solve(crossprod(mZ,tryCatch(solve(mD), error = function(e) ginv(mD)))%*%mZ + diag(k)),  error = function(e) ginv(crossprod(mZ,tryCatch(solve(mD), error = function(e) ginv(mD)))%*%mZ + diag(k)))
    mM <- tcrossprod(mS,as.matrix(mZ))%*%tcrossprod(tryCatch(solve(mD), error = function(e) ginv(mD)),y-mMu) # note, we have a matrix of m, each column for each observation
    
    eU <- mM
    
    eUU.sum <- matrix(0, nrow=k, ncol=k)
    for(j in 1:n){
      eUU.sum <-eUU.sum +  mS +tcrossprod(eU[,j],eU[,j])
    }
    
    eU.sum <- matrix(0,nrow=k, ncol=p)
    for(j in 1:n){
      eU.sum <-eU.sum+ tcrossprod(eU[,j], y[j,]-vmu)  
      
    }
    
    mZ.new <- crossprod(eU.sum, tryCatch(solve(eUU.sum), error = function(e) ginv(eUU.sum)))
    mD.new <- diag(diag(cov(y)*((n-1)/n) - (mZ.new%*%eU.sum/n))) #note I have corrected the sample cov matrix to data cov matrix
    
    crit <- mean((mD - mD.new)^2)
    
    mD <- mD.new
    mZ <- mZ.new
    
    
  }
  
  
  return(list(mD=mD, mZ=mZ, mM=mM))
}


########################################################################################

## Function to initialise genDA starting values ##

.genDA_start_values<- function(y, X = NULL, family, num.lv = 2) {
  
  set.seed(1)

  N <-n <- nrow(y); 
  m <- ncol(y); 
  y <- as.matrix(y)
  
  num.X <- 0; 
  if(!is.null(X)) num.X <- dim(X)[2]
  
  if(!is.numeric(y))
    stop("y must be numeric")
  
  if(num.lv > 0) {
    mU <- mvtnorm::rmvnorm(N, rep(0, num.lv));
  }
  
  if(num.lv == 0) { mU <- NULL }

  if(length(family)==1 & m > 1){
    
    if(!(family %in% c("poisson","binomial", "gaussian", "log-normal", "negative-binomial", "ZIP")))
      stop("Input family not supported")
    
    if(family=="log-normal"){
      family="gaussian"
      y <- log(y)
    }
    
    if(family=="negative-binomial"){
      family="negative.binomial"
    }
    
    if(family =="ZIP"){
      family = "poisson"
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
      
      if(response_types[j]=="ZIP"){
        response_types[j]="poisson"
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
  if(num.lv > 0){
    if(m > 2 && n > 2){
      
      fa <-  .EM.FA(resi, num.lv)
      mLambda <- (fa$mZ)
      mU <- t(fa$mM)
      
    } else {
      mLambda <- matrix(1,m,num.lv)
      mLambda[upper.tri(mLambda)]=0
      mU <- matrix(0,n,num.lv)
    }
  }
  params <- mLambda
  
  if(num.lv > 0) {
    mU <- mU+mvtnorm::rmvnorm(n, rep(0, num.lv),diag(num.lv))
  }

  out <- list(mLambda=mLambda)
  
  if(num.lv > 0) out$mU <- mU
  
  return(out)
}

#######################################################################################

.vec2mat <- function(vy) {
  vy=as.numeric(vy)
  K <- length(unique(vy))
  n <- length(vy)
  mY <- matrix(0, nrow = n, ncol = K)
  for (i in 1:length(vy)) {
    val <- vy[i]
    mY[i, val] <- 1
  }
  return(mY)
}

########################################################################################


.mat2vec <- function(mY) {
  classval <- c()
  for (i in 1:nrow(mY)) {
    v <- mY[i, ]
    classval[i] <- which(v == 1)
  }
  return(classval)
}

#######################################################################################

quiet <- function(x) { 
       sink(tempfile()) 
       on.exit(sink()) 
       invisible(force(x)) 
}

#######################################################################################

.logMatToGamma <- function(log.mA) {
  R <- nrow(log.mA)
  C <- ncol(log.mA)
  vm <- matrixStats::rowMaxs(log.mA)
  mM <- matrix(vm, R, C, byrow = FALSE)
  mA.til <- exp(log.mA - mM)
  vs <- matrixStats::rowSums2(mA.til)
  mS <- matrix(vs, R, C, byrow = FALSE)
  mGamma <- mA.til / mS
  return(mGamma)
}

#######################################################################################

.num.2.fac <- function(x, vy.fac){
  dat <- as.numeric(unique(vy.fac))
  names(dat) <- unique(vy.fac)
  val=names(which(dat==x))
  return(val)
}

