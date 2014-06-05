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
VectorXd FunctionSetup::vector_u(){
  VectorXd u(this->A.cols()); 
  u= this->A.transpose()*this->mean;
  return u;
}

/// Generates C(Sigma) = A^T*(beta*betaI + SIGMA)*A

MatrixXd FunctionSetup::matrix_C(){
  double betaSquare = this->beta*this->beta;
  MatrixXd C(this->A.rows(),this->A.cols());
  this->covar.diagonal().array() + betaSquare;
  C = A.transpose() * this->covar; 
  C = C * A;
  return C;
}

///Returns value of A*C(SIGMA)^(-1)*z - A*C(SIGMA)^(-1)*u(mu)
VectorXd FunctionSetup::function_g(VectorXd z){
  VectorXd result(z.size());
  result = (A*this->matrix_C().inverse()*z - A * this->matrix_C().inverse()*this->vector_u());
  return result; 
}

int main(void){
  
}