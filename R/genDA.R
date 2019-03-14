library(TMB)

setwd("C:/Users/sarah/Dropbox (Sydney Uni)/Uni/PhD/Projects/genDA-chapter/Scripts/DoubleCode")
compile("genDA_f.cpp")
dyn.load(dynlib("genDA_f"))

### Simulate Data ###

# Generate data

n <- 100
m <- 10
d <- 2
p <- 1

vz <- rbinom(n,1,0.5)
mX_cov <- matrix(vz,n,p)

vtau <- rep(0.0,n)
vbeta0 <- rep(-1,m)
mB <- matrix(rep(2,m),p,m)
response_types <- rep(c("logistic"),each=m)

mLambda <- matrix(rnorm(m*d,sd=0.1),m,d) 
mLambda[upper.tri(mLambda)] <- 0
for (k in 1:d) {
  mLambda[k,k] <- 1
}
mL <- t(mLambda)

mU <- matrix(0,n,d)
for (i in 1:n) {
  vu <- rnorm(d,sd=1)
  mU[i,] <- vu
}

mEta.fixed  <- matrix(vtau)%*%matrix(1,1,m) + matrix(1,n,1)%*%vbeta0 + mX_cov%*%mB
mEta.latent <- mU%*%mL

mEta <- mEta.fixed + mEta.latent


mY <- matrix(0,n,m)
for (j in 1:m) {
  if (response_types[j]=="logistic") {
    mY[,j] <- rbinom(n,1,1/(1+exp(-mEta[,j])))
   
  }
  if (response_types[j]=="Poisson") {
    mY[,j] <- rpois(n,exp(mEta[,j]))
    
  }
}


##############################################

# How do we do if we just use glms

vbeta0.hat <- c()
mB.hat <- matrix(0,p,m)
for (j in 1:m) {
  
  if (response_types[j]=="logistic") {
    res <- glm(mY[,j]~mX_cov,family=binomial)
    response_types[j]=1
  }
  if (response_types[j]=="Poisson") {
    res <- glm(mY[,j]~mX_cov,family=poisson)
    response_types[j]=2
  }
  
  vbeta0.hat[j] <- res$coef[1]
  mB.hat[,j] <- res$coef[-1]
}

response_types <- as.numeric(response_types)

################################################################################

# Set Ridge regression constants

sigma2_beta0	<- 1.0E0
vsigma2_beta    <- rep(1.0E0,m)
vsigma2_lambda  <- rep(1.0E2,m)
vsigma2_tau     <- rep(1.0E-2,n)

vsigma2 <- c()
vsigma2[1]           <- sigma2_beta0
vsigma2[1+(1:m)]     <- vsigma2_beta
vsigma2[1+m+(1:m)]   <- vsigma2_lambda
vsigma2[1+m+m+(1:n)] <- vsigma2_tau

vtau.init <- 0*vtau
vbeta0.init <- vbeta0.hat
mB.init <- mB.hat


mL.init <- 0*(matrix(0.01*rnorm(m*d),d,m) + mL)
for (k in 1:d) {
  mL.init[k,k] <- 1
}

mU.init <- matrix(rnorm(n*d),n,d) 

vtheta <- c(vtau.init, vbeta0.init, mB.init, mU.init, mL.init)

data <- list()
data$mY <- mY
data$mX_cov <- mX_cov
data$vsigma2 <- vsigma2
data$response_types <- response_types
data$d <- d

parameters <- list()
parameters$vtheta <- vtheta

obj = MakeADFun(data,parameters,DLL = "genDA_f")
#opt = optim(obj$par,obj$fn, obj$gr)

vals <-  obj$report()

vals$sigma2_beta0
sigma2_beta0

head(vals$vsigma2_beta)
head(vsigma2_beta)

head(vals$vsigma2_lambda)
head(vsigma2_lambda)

head(vals$vsigma_tau)
head(vsigma2_tau)

tail(vals$vtau)
tail(vtau)

tail(vals$vbeta0)
tail(vbeta0.init)

vals$mB
mB.init

head(vals$mU)
head(mU.init)
