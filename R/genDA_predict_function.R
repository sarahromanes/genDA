genDA_predict <- function(mY_test, res_genDA, prior_beta = c(1,1)){
  
  
  mB <- res_genDA$mB
  if(is.null(mB)){
    stop("GLLVM not trained with class responses included. Try fitting genDA again with class information captured in mX.")
  }
  
  n_test <- nrow(mY_test)
  p <- res_genDA$p
  m <- ncol(mY_test)
  d <- res_genDA$d
  family <- res_genDA$family
  n_train <- res_genDA$n
  
  A <- prior_beta[1]
  B <- prior_beta[2]

  vsigma2_tau     <- rep(1.0E-4,n_test)
  vtau.init <- rep(0,n_test)
  
  mU.init <- matrix(rnorm(n_test*d), nrow=n_test, ncol=d) # get same results if I use start.values.gllvm.TMB or random. 
  
  mX_pred <- rep(0, n_test)
  prob_X <- matrix(0,nrow = n_test, ncol=2)
  colnames(prob_X) <- c("Class 0", "Class 1")
  
  for(i in 1:n_test){
    
    data_1 <- list()
    data_1$mY <- t(as.matrix(mY_test[i, ]))
    data_1$mX <- t(as.matrix(c(1), nrow=1, ncol=p)) #setting vc_test[i] = 1
    data_1$vsigma2_tau <- vsigma2_tau[i]
    data_1$response_types <- res_genDA$tmb_types
    data_1$d <- d
    data_1$mB <- res_genDA$mB
    data_1$mL <- res_genDA$mL
    data_1$vbeta0 <- as.matrix(res_genDA$vbeta0)
    data_1$vphi <- as.matrix(res_genDA$vphi)
    
    data_0 <- data_1
    data_0$mX <-  1- data_1$mX # same parameters as before, however now vc_test[i]=0
        
    parameters <- list()
    parameters$mU <- t(as.matrix(mU.init[i,]))
    parameters$vtau <- vtau.init[i]
    
    
    obj_1 <- MakeADFun(data_1,parameters,DLL = "genDA_f_predict_approximation", silent=TRUE)
    obj_0 <-  MakeADFun(data_0,parameters,DLL = "genDA_f_predict_approximation", silent=TRUE)
    opt_1 <-  nlminb(obj_1$par,obj_1$fn, obj_1$gr,control=list(rel.tol=1.0E-6))
    opt_0 <- nlminb(obj_0$par,obj_0$fn, obj_0$gr,control=list(rel.tol=1.0E-6))
    
    ll_1 <- -1*(opt_1$objective) +  lbeta(A + sum(res_genDA$vc) + 1, B + n_train - sum(res_genDA$vc) -1) 
    ll_0 <- -1*(opt_0$objective) +  lbeta(A + sum(res_genDA$vc), B + n_train - sum(res_genDA$vc)) 
    
    if(ll_1 > ll_0 ){
      mX_pred[i] = 1
    }
    
    prob_X[i, 1] <- exp(ll_0)/(exp(ll_0)+exp(ll_1))
    prob_X[i, 2] <- 1- prob_X[i,1]
    
  }
  
  return(list(mX_pred = mX_pred, 
              prob_X = round(prob_X,4)))
  
}