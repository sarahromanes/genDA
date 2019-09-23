// --- genDA_main.cpp ---

#define TMB_LIB_INIT R_init_genDA

#include <TMB.hpp>
#include "genDA_f.hpp"
#include "genDA_f_null_dis.hpp"
#include "genDA_f_null_X.hpp"
#include "genDA_f_null_X_dis.hpp"
#include "genDA_f_predict.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model_name);
  if(model_name == "genDA_f") {
    return genDA_f(this);
  } else if(model_name == "genDA_f_null_dis") {
    return genDA_f_null_dis(this);
  } else if(model_name == "genDA_f_null_X"){
      return genDA_f_null_X(this);
  } else if(model_name == "genDA_f_null_X_dis"){
      return genDA_f_null_X_dis(this);
  } else if(model_name =="genDA_f_predict"){
    return genDA_f_predict(this);
  } else {
    error("Unknown model_name.");
  }
  return 0;
}
