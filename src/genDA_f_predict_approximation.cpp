#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(mY);
  DATA_MATRIX(mX);
  DATA_SCALAR(vsigma2_tau);
  DATA_IVECTOR(response_types); // response_types needs to be coded in as integer 1 (Bernoulli) or 2 (Poisson), or 3 (Gaussian), or 4 (Log-Normal)
  DATA_INTEGER(d);
  DATA_MATRIX(mB); // optimised from previous TMB call
  DATA_MATRIX(mL); // optimised from previous TMB call
  DATA_MATRIX(vbeta0); //optimised from previous TMB call
  DATA_MATRIX(vsigma_norm); //optimised from previous TMB call
  PARAMETER_MATRIX(mU);
  PARAMETER(vtau);

  int n = 1 ;
  int m = mY.array().cols();

  matrix<Type> mOnes(n,m); // create a matrix of ones, and then extract a row and column below to get mOnes (length m, and n) below
  for(int i =0; i < n; i++){
      for(int j = 0; j < m; j++){
          mOnes(i,j) = 1;
      }
  }

  matrix<Type> oneM = mOnes.array().row(0);
  

 matrix<Type> mEta = vtau*oneM + vbeta0.transpose() + mX*mB +  mU*mL;

  REPORT(oneM);
  REPORT(mX);
  REPORT(mB);
  REPORT(mU);
  REPORT(mY);
  REPORT(vtau);
  REPORT(vbeta0);


    // CALCULATE LOG LIKELIHOOD

   Type nll = 0; // (We will use  += -= *= /= incremental operators to add/subtract to the nll in increments. Note this is the negative ll and as such we will subtract when we mean to add etc... )

    for(int j = 0; j < m; j++){
      if(response_types(j)==1){ 
        // BERNOULLI DISTRIBUTION
        nll -= mY(0,j)*mEta(0,j) - log(1 + exp(mEta(0,j)));
      }
      if(response_types(j)==2){
        // POISSON DISTRIBUTION
        nll-=  dpois(mY(0,j), exp(mEta(0,j)), true);
      }
      if(response_types(j)==3){
        // GAUSSIAN DISTRIBUTION
        nll -= dnorm(mY(0,j), mEta(0,j), vsigma_norm(j), true);
      }
      if(response_types(j)==4){
        // LOG NORMAL DISTRIBUTION
        nll -= dnorm(log(mY(0,j)), mEta(0,j), vsigma_norm(j), true);
      }
    }
  
  REPORT(mEta);
  
  // // CALCULATE AND "ADD" REGULARISATION TERMS 
  
  nll += 0.5*pow(vtau,2.0)/vsigma2_tau;
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
        Type sderiv = (exp(mEta(i,j))*((1.0 + exp(mEta(i,j))) -1.0))/( pow((1.0 + exp(mEta(i,j))), 2.0));
        mVal += sderiv*(lambda_j*lambda_j.transpose());
      }
      if(response_types(j)==2){
        // POISSON DISTRIBUTION
        mVal += exp(mEta(i,j))*(lambda_j*lambda_j.transpose());
      }
      if(response_types(j)==3){
        // GAUSSIAN DISTRIBUTION
        mVal += (1/(pow(vsigma_norm(j),2.0)))*(lambda_j*lambda_j.transpose());
      }
      if(response_types(j)==4){
        // LOGNORMAL DISTRIBUTION
        mVal += (1/(pow(vsigma_norm(j),2.0)))*(lambda_j*lambda_j.transpose());
      }
    }
    matrix<Type> mDet = mVal + mI;
    
    nll += 0.5*log(mDet.determinant());
  }
  
  return(nll);
 
}
