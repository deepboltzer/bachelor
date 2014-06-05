#include <eigen3/Eigen/Dense>
#include "FunctionSetup.h"

using namespace Eigen;

/// Default Constructor.
FunctionSetup::FunctionSetup(){
  
}

/// Constructor which generates FunctionSetup with param mean,CovarianceMatrix, Matrix A and trueskill factor beta.
FunctionSetup::FunctionSetup(VectorXd m, MatrixXd S, MatrixXd partial, double b)
: mean(m), covar(S), A(partial), beta(b)
{
  
}

/// Generates u(mu) = A^T*mu as VectorXd.
VectorXd FunctionSetup::function_u(){
  VectorXd u(this->A.cols()); 
  u= this->A.transpose()*this->mean;
  return u;
}

int main(void){
  
}