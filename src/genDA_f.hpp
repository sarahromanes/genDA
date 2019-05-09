// --- genDA_f.hpp ---

// genDA_f is the general genDA algorithm used to fit GLLVMs when there are covariates to be estimated.

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template <class Type>
Type genDA_f(objective_function<Type>* obj) {

  DATA_MATRIX(y);
  DATA_MATRIX(X);
  DATA_SCALAR(sigma2_beta0);
  DATA_VECTOR(vsigma2_beta);
  DATA_VECTOR(vsigma2_lambda);
  DATA_VECTOR(vsigma2_tau);
  DATA_IVECTOR(response_types); // response_types needs to be coded in as integer 1 (Bernoulli) or 2 (Poisson), or 3 (Gaussian), or 4 (Log-Normal), or 5 (NB)
  DATA_INTEGER(d);
  DATA_IVECTOR(vphi_inds);
  
  PARAMETER_VECTOR(log_vphi);
  PARAMETER_MATRIX(mB);
  PARAMETER_VECTOR(lambda);
  PARAMETER_MATRIX(mU);
  PARAMETER_MATRIX(vtau);
  PARAMETER_MATRIX(vbeta0);

  int n = y.array().rows();
  int m = y.array().cols();

  vector<Type> vphi = exp(log_vphi);
    
  // To create lambda as matrix upper triangle

  matrix<Type> mL(d,m);
  if(d>0){

    for (int j=0; j<m; j++){
      for (int k=0; k<d; k++){
        if (j < k){
          mL(k,j) = 0;
        } else{
          mL(k,j) = lambda(j);
          if (k > 0){
            mL(k,j) = lambda(k+j+k*m-(k*(k-1))/2-2*k);
          }
        }
      }
    }
  }
  
  matrix<Type> mOnes(n,m); // create a matrix of ones, and then extract a row and column below to get mOnes (length m, and n) below
  for(int i =0; i < n; i++){
      for(int j = 0; j < m; j++){
          mOnes(i,j) = 1;
      }
  }

  matrix<Type> oneM = mOnes.array().row(0);
  matrix<Type> oneN = mOnes.array().col(0);

  matrix<Type> mEta = vtau*oneM + oneN*vbeta0.transpose() + X*mB +  mU*mL;

    // CALCULATE LOG LIKELIHOOD

  Type nll = 0; // (We will use  += -= *= /= incremental operators to add/subtract to the nll in increments. Note this is the negative ll and as such we will subtract when we mean to add etc... )

  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      if(response_types(j)==1){ 
        // BERNOULLI DISTRIBUTION
        nll -= y(i,j)*mEta(i,j) - log(1 + exp(mEta(i,j)));
      }
      if(response_types(j)==2){
        // POISSON DISTRIBUTION
        nll-=  dpois(y(i,j), exp(mEta(i,j)), true);
      }
      if(response_types(j)==3){
        // GAUSSIAN DISTRIBUTION
        nll -= dnorm(y(i,j), mEta(i,j), vphi(vphi_inds(j)), true);
      }
      if(response_types(j)==4){
        // LOG NORMAL DISTRIBUTION
        nll-= dnorm(log(y(i,j)), mEta(i,j), vphi(vphi_inds(j)), true);
      }
      if(response_types(j)==5){
        // NEGATIVE BINOMIAL DISTRIBUTION
        Type mu = exp(mEta(i,j));
        Type var = mu + pow(mu,2.0)/vphi(vphi_inds(j));
        nll-= dnbinom2(y(i,j), mu, var, true);
      }
    }
  }

  // CALCULATE AND "ADD" REGULARISATION TERMS 
  
  nll += 0.5*(vtau.array().pow(2.0)/vsigma2_tau.array()).array().sum(); 
  nll += 0.5*(vbeta0.array().pow(2.0)/sigma2_beta0).array().sum();
  
  for(int j = 0; j < m; j++){
    nll += 0.5*(mB.array().col(j).pow(2.0)/vsigma2_beta(j)).array().sum();
    nll += 0.5*(mL.array().col(j).pow(2.0)/vsigma2_lambda(j)).array().sum();
  }
  
  for(int i = 0; i < n; i++){
    nll += 0.5*(mU.array().row(i).pow(2.0)).array().sum();
  }
  
  
  // CALCULATE AND "ADD" LOG- DETERMINANT TERM

  // Create d x d identity matrix

  matrix<Type> mI(d,d);
  mI.setZero(); 
  
  for(int i = 0; i < d; i++){
    for(int j = 0; j < d; j++){
      if(i == j){
        mI(i,j)=1.0;
      }
    }
  }

  // Then iterate through i and to calculate log-determinant term

  for(int i = 0; i < n; i++){
    matrix<Type> mVal(d,d);
    mVal.setZero(); 
    for(int j = 0; j < m; j++){
      matrix<Type> lambda_j = mL.array().col(j);
      if(response_types(j)==1){
        // BERNOULLI DISTRIBUTION
        Type sderiv = (exp(mEta(i,j))*((1.0 + exp(mEta(i,j))) -1.0))/( pow((1.0 + exp(mEta(i,j))), 2.0));
        mVal += sderiv*(lambda_j*lambda_j.transpose());
      }
      if(response_types(j)==2){
        // POISSON DISTRIBUTION
        mVal += exp(mEta(i,j))*(lambda_j*lambda_j.transpose());
      }
      if(response_types(j)==3){
        // GAUSSIAN DISTRIBUTION
        mVal += (1/(pow(vphi(vphi_inds(j)),2.0)))*(lambda_j*lambda_j.transpose());
      }
      if(response_types(j)==4){
        // LOGNORMAL DISTRIBUTION
        mVal += (1/(pow(vphi(vphi_inds(j)),2.0)))*(lambda_j*lambda_j.transpose());
      }
      if(response_types(j)==5){
        // NEGATIVE BINOMIAL DISTRIBUTION
        Type sderiv = ((y(i,j)+vphi(vphi_inds(j)))*(vphi(vphi_inds(j))*exp(mEta(i,j)))/pow(vphi(vphi_inds(j)) + exp(mEta(i,j)),2.0));
        mVal += sderiv*(lambda_j*lambda_j.transpose());
      }
    }
    matrix<Type> mDet = mVal + mI;
    
    nll += 0.5*log(mDet.determinant());
  }
  
  return(nll);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
