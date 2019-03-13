#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(mY);
  DATA_MATRIX(mX_cov);
  DATA_MATRIX(mOnes);
  DATA_VECTOR(vsigma2);
  DATA_IVECTOR(respose_types); // response_types needs to be coded in as integer 1 (Bernoulli) or 2 (Poisson)
  DATA_INTEGER(d);
  PARAMETER_VECTOR(vtheta);
    
  int n = mX_cov.array().rows();
  int m = mY.array().cols();
  int p = mX_cov.array().cols();
  
  Type sigma2_beta0 = vsigma2(0); // remembering that C++ starts at index 0
  matrix<Type> vsigma2_beta = vsigma2.segment(1,m);
  matrix<Type> vsigma2_lambda = vsigma2.segment(m+1, m+m);
  matrix<Type> vsigma_tau = vsigma2.segment(m+m+1, n+m+m);
  
  matrix<Type> vtau = vtheta.segment(0, n-1);
  matrix<Type> vbeta0 = vtheta.segment(n, n+m-1);
  
  matrix<Type> mB(p,m); 
  for(int j = 0; j < m; j++){
    mB.array().col(j) = vtheta.segment(n+m+((j-1)*p), p+ n+m+((j-1)*p)-1); // update the jth col of mB
  }
  
  matrix<Type> mU(n,d);
  for(int i = 0; i <n; i++){
    mU.array().row(i)= vtheta.segment(n+m+((i-1)*d), d + n+m+((i-1)*d) -1); // update the ith row of mU
  }
  
  matrix<Type> mL(d,m);
  mL.setZero();  // Since vectors, matrices and arrays are not zero-initialized in C++, a zero initialized object is created using Eigens setZero():
  
  int count = n + m + (m*p) +(n*d);
  for(int j = 0; j < m; j++){ // two layers of iteration are needed - was not possible to subset cols/rows of matrices in C++
    for(int k = 0; k < d; k++){
      if(j >= k){
        mL(j,k) = vtheta(count + (d*j) + k);
      }
      if(j == k){
        mL(j, k) = exp(mL(j,k));
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
  
  matrix<Type> mEta_fixed = vtau.transpose()*oneM + oneN.transpose()*vbeta0 + mX_cov*mB;
  matrix<Type> mEta_latent = mU*mL;
  matrix<Type> mEta = mEta_fixed + mEta_latent;
  
  // CALCULATE LOG LIKELIHOOD
  matrix<Type> ml(n,m);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      if(response_types(j)==1){ 
        // BERNOULLI DISTRIBUTION
        ml(i,j) = mY(i,j)*mEta(i,j) - log(1 + exp(mEta(i,j)));
      }
      if(response_types(j)==2){
        // POISSON DISTRIBUTION
        ml(i,j) = mY(i,j)*mEta(i,j) -  exp(mEta(i,j));
      }
    }
  }
  
  matrix<Type> vl = ml.colwise().sum();
  
  // CALCULATE REGULARISATION TERMS 
  Type pen = 0;
  
  pen -= 0.5*(vtau.array().pow(2.0)/vsigma_tau.array()).array().sum(); //  += -= *= /= incremental operators
  pen -= 0.5*(vbeta0.array().pow(2.0)/sigma2_beta0).array().sum();
  
  for(int j = 0; j < m; j++){
    pen -= 0.5*(mB.array().col(j).pow(2.0)/vsigma2_beta(j)).array().sum();
    pen -= 0.5*(mL.array().col(j).pow(2.0)/vsigma2_lambda(j)).array().sum();
  }
  
  for(int i = 0; i < n; i++){
    pen -= 0.5*(mU.array().row(i).pow(2.0)).array().sum();
  }

  // COMBINE INTO OVERALL PENALISED NLL

  Type nll = -(sum(vl) + pen); //return negative of total value (not sure if you can fn scale this? just in case - isn't too hard anyways)

  return(nll);
}
