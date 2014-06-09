#include "../gaussians/gaussian.h"
#include <math.h>
#include "../gaussians/constants.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Cholesky>
#include <stdlib.h>
#include <iostream>

using namespace Eigen;
/// TODO: INCLUDE FUNCTION g from function setup:

/// Returns aprroximation of an truncated gaussian with mean mu and covariance covar. 
VectorXd ApproximationAlgorithm(VectorXd mu, MatrixXd covar, VectorXd lower, VectorXd upper, int N, VectorXd w ){
  VectorXd approx(mu.size());
  LLT<MatrixXd> lltOfA(covar); // compute the Cholesky decomposition of covar
  MatrixXd L = lltOfA.matrixL(); // retrieve factor L  in the decomposition
  int n = mu.size();
  VectorXd z(n);
  double u,l;
  Gaussian G1(0,1);
  for(int k = 0;k<N;k++){
    for(int i =0;i<n;i++){
      double nominator,denominator,sum;
      for(int j = 0;j<i-1;j++){
	sum += L(i,j)*z[j];
      }
      nominator = lower[i] - mu[i] -sum;
      denominator = L(i,i);
      l = G1.cumulativeTo(nominator/denominator);
      nominator = upper[i] - mu[i] -sum;
      u = G1.cumulativeTo(nominator/denominator);
      sum = 0;
      z[i] = G1.inverseCumulativeTo(l+ w[i]*(u-l));
    }
    approx = (k/(k+1)) * approx + (1/(k+1)) * (L*z + mu); 
  }
  return approx;
}

int main(void){
  
}