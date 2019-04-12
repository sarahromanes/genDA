EM.FA <- function(mY, k){
  
  n <- nrow(mY)
  p <- ncol(mY)
  
  pc <- prcomp(mY)
  mD <- diag(p)
  mZ <- pc$rotation[,1:k]
  
  
  # Our estimate of vmu is calculated once out the iteration loop.
  
  vmu <- matrixStats::colMeans2(mY)
  mMu <- matrix(rep(vmu, n), nrow=n, ncol=p, byrow=T)
  
  min.err <- 1.0E-8
  
  crit <- 1
  
  
  while(crit > min.err){
    
    mS <- solve(crossprod(mZ,solve(mD))%*%mZ + diag(k))
    mM <- tcrossprod(mS,as.matrix(mZ))%*%tcrossprod(solve(mD),mY-mMu) # note, we have a matrix of m, each column for each observation
    
    eU <- mM
    
    eUU.sum <- matrix(0, nrow=k, ncol=k)
    for(j in 1:n){
      eUU.sum <-eUU.sum +  mS +tcrossprod(eU[,j],eU[,j])
    }
    
    eU.sum <- matrix(0,nrow=k, ncol=p)
    for(j in 1:n){
      eU.sum <-eU.sum+ tcrossprod(eU[,j], mY[j,]-vmu)  
      
    }
    
    mZ.new <- crossprod(eU.sum,solve(eUU.sum))
    mD.new <- diag(diag(cov(mY)*((n-1)/n) - (mZ.new%*%eU.sum/n))) #note I have corrected the sample cov matrix to data cov matrix
    
    crit <- mean((mD - mD.new)^2)
    
    mD <- mD.new
    mZ <- mZ.new
    
    
  }
  
  
  return(list(mD=mD, mZ=mZ, mM=mM))
}