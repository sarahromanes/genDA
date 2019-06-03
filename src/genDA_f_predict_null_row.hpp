// --- genDA_f_predict.hpp ---

// genDA_f_predict is the objective function for the predictive component of genDA. It requires that a GLLVM has been fit using class information

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template <class Type>
Type genDA_f_predict_null_row(objective_function<Type>* obj) {

  DATA_MATRIX(y);
  DATA_MATRIX(X);
  DATA_IVECTOR(response_types); 
  DATA_INTEGER(d);
  DATA_MATRIX(mB); // optimised from previous TMB call
  DATA_MATRIX(mL); // optimised from previous TMB call
  DATA_MATRIX(vbeta0); //optimised from previous TMB call
  DATA_MATRIX(vphi); //optimised from previous TMB call
  PARAMETER_MATRIX(mU);

  int n = 1 ;
  int m = y.array().cols();

  matrix<Type> mOnes(n,m); // create a matrix of ones, and then extract a row and column below to get mOnes (length m, and n) below
  for(int i =0; i < n; i++){
      for(int j = 0; j < m; j++){
          mOnes(i,j) = 1;
      }
  }

  matrix<Type> oneM = mOnes.array().row(0);
  
  matrix<Type> mEta = vbeta0.transpose() + X*mB +  mU*mL;


    // CALCULATE LOG LIKELIHOOD

   Type nll = 0; // (We will use  += -= *= /= incremental operators to add/subtract to the nll in increments. Note this is the negative ll and as such we will subtract when we mean to add etc... )

    for(int j = 0; j < m; j++){
      if(response_types(j)==1){ 
        // BERNOULLI DISTRIBUTION
        nll -= y(0,j)*mEta(0,j) - log(1 + exp(mEta(0,j)));
      }
      if(response_types(j)==2){
        // POISSON DISTRIBUTION
        nll-=  dpois(y(0,j), exp(mEta(0,j)), true);
      }
      if(response_types(j)==3){
        // GAUSSIAN DISTRIBUTION
        nll -= dnorm(y(0,j), mEta(0,j), vphi(j), true);
      }
      if(response_types(j)==4){
        // LOG NORMAL DISTRIBUTION
        nll -= dnorm(log(y(0,j)), mEta(0,j), vphi(j), true);
      }
      if(response_types(j)==5){
        // NEGATIVE BINOMIAL DISTRIBUTION
        Type mu = exp(mEta(0,j));
        Type var = mu + pow(mu,2.0)/vphi(j);
        nll-= dnbinom2(y(0,j), mu, var, true);
      }
      if(response_types(j)==6){
        // ZERO INFLATED POISSON DISTRIBUTION
        vphi(j) = vphi(j)/ (1.0 +vphi(j));
        nll -= dzipois(y(0,j), exp(mEta(0,j)),vphi(j), true); 
      }
      
    }
  
  REPORT(mEta);
  
  // // CALCULATE AND "ADD" REGULARISATION TERMS 
  
  nll += 0.5*(mU.array().row(0).pow(2.0)).array().sum();
  
  
  // // // CALCULATE AND "ADD" LOG- DETERMINANT TERM

  // // // Create d x d identity matrix

  matrix<Type> mI(d,d);
  mI.setZero(); 
  
  for(int i = 0; i < d; i++){
    for(int j = 0; j < d; j++){
      if(i == j){
        mI(i,j)=1.0;
      }
    }
  }

  // // // // Then iterate through i and to calculate log-determinant term

  for(int i = 0; i < n; i++){
    matrix<Type> mVal(d,d);
    mVal.setZero(); 
    for(int j = 0; j < m; j++){
      matrix<Type> lambda_j = mL.array().col(j);
      if(response_types(j)==1){
        // BERNOULLI DISTRIBUTION
        Type sderiv = 1.0 /( pow((1.0 + exp(mEta(i,j))), 2.0));
        mVal += sderiv*(lambda_j*lambda_j.transpose());
      }
      if(response_types(j)==2){
        // POISSON DISTRIBUTION
        mVal += exp(mEta(i,j))*(lambda_j*lambda_j.transpose());
      }
      if(response_types(j)==3){
        // GAUSSIAN DISTRIBUTION
        mVal += (1/(pow(vphi(j),2.0)))*(lambda_j*lambda_j.transpose());
      }
      if(response_types(j)==4){
        // LOGNORMAL DISTRIBUTION
        mVal += (1/(pow(vphi(j),2.0)))*(lambda_j*lambda_j.transpose());
      }
      if(response_types(j)==5){
        // NEGATIVE BINOMIAL DISTRIBUTION
        Type sderiv = ((y(i,j)+vphi(j))*(vphi(j)*exp(mEta(i,j)))/pow(vphi(j) + exp(mEta(i,j)),2.0));
        mVal += sderiv*(lambda_j*lambda_j.transpose());
      }
        if(response_types(j)==6){
        // ZERO INFLATED POISSON DISTRIBUTION
        Type sderiv = 0.0;
        if(y(i,j)==0){ 
         sderiv = -1*( (vphi(j) -1.0)*exp(mEta(i,j))*(vphi(j)*(-exp(exp(mEta(i,j))) + exp(mEta(i,j) + exp(mEta(i,j))) + 1 ) -1 ))/pow(( (vphi(j)*(exp(exp(mEta(i,j))) -1.0 ) ) + 1.0 ), 2.0);
        } else {
        sderiv = exp(mEta(i,j));
        }
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
