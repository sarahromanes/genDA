#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(mY);
  DATA_MATRIX(mX);
  DATA_VECTOR(vsigma2);
  DATA_IVECTOR(response_types); // response_types needs to be coded in as integer 1 (Bernoulli) or 2 (Poisson), or 3 (Gaussian), or 4 (Log-Normal)
  DATA_INTEGER(d);
  PARAMETER_VECTOR(log_vsigma_norm);
  PARAMETER_MATRIX(mB);
  PARAMETER_VECTOR(lambda);
  PARAMETER_MATRIX(mU);
  PARAMETER_MATRIX(vtau);
  PARAMETER_MATRIX(vbeta0);

  int n = mX.array().rows();
  int m = mY.array().cols();

  vector<Type> vsigma_norm = exp(log_vsigma_norm);
  
  Type sigma2_beta0 = vsigma2(0); // remembering that C++ starts at index 0
  matrix<Type> vsigma2_beta = vsigma2.segment(1,m); // vector.segement(i,n) takes a Block containing n elements, starting at position i.
  matrix<Type> vsigma2_lambda = vsigma2.segment(m+1, m); 
  matrix<Type> vsigma_tau = vsigma2.tail(n); // vector.tail(n) takes a Block containing the last n elements
  
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

  matrix<Type> mEta = vtau*oneM + oneN*vbeta0.transpose() + mX*mB +  mU*mL;

    // CALCULATE LOG LIKELIHOOD

  Type nll = 0; // (We will use  += -= *= /= incremental operators to add/subtract to the nll in increments. Note this is the negative ll and as such we will subtract when we mean to add etc... )

  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      if(response_types(j)==1){ 
        // BERNOULLI DISTRIBUTION
        nll -= mY(i,j)*mEta(i,j) - log(1 + exp(mEta(i,j)));
      }
      if(response_types(j)==2){
        // POISSON DISTRIBUTION
        nll-=  dpois(mY(i,j), exp(mEta(i,j)), true);
      }
      if(response_types(j)==3){
        // GAUSSIAN DISTRIBUTION
        nll -= dnorm(mY(i,j), mEta(i,j), vsigma_norm(j), true);
      }
      if(response_types(j)==4){
        // LOG NORMAL DISTRIBUTION
        nll-= dnorm(log(mY(i,j)), mEta(i,j), vsigma_norm(j), true);
      }
    }
  }

  // CALCULATE AND "ADD" REGULARISATION TERMS 
  
  nll += 0.5*(vtau.array().pow(2.0)/vsigma_tau.array()).array().sum(); 
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
  

  // REPORT VALUES TO R FOR CHECKING

  REPORT(mL);
  REPORT(mU);
  REPORT(mB);
  REPORT(vtau);
  REPORT(vbeta0);
  REPORT(sigma2_beta0);
  REPORT(vsigma2_lambda);
  REPORT(vsigma2_beta);
  REPORT(vsigma_tau);
  REPORT(mEta);
  REPORT(nll);

  return(nll);
 
}
