#include <math.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Cholesky>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "gameSetup.h"
#include "../gaussians/gaussian.h"

/// TODO: INCLUDE FUNCTION g from function setup:
/// At the moment we assume that g is constant function, as in Ranking Probability Algorithm.
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

/// TODO: Implement rest of update Algorithm. Need more information about function g. 
void updateAlgortithm(gameSetup S, double drawMargin){
  /// Set number players and number teams:
  int n = S.getNumberPlayers();
  int k = S.getNumberTeams();
  /// Integration limit vectors:
  vector<double> a;
  vector<double> b;
  /// First Update of player variance:
  for(int j = 0; j< k; j++){
    double value = sqrt(S.getSigmaByIndex(j)*S.getSigmaByIndex(j) + S.getTau()*S.getTau());
    S.setSigma(j,value);
  }
  
  /// Declaration of PlayerTeamMatrix A: 
  MatrixXd A(n,k-1);
  
  /// Set values of matrix A:
  bool first = true;
  int lower1 = 0;
  int upper1 = 0;
  int lower2 = 0;
  int upper2 = 0;
  
  for(int j = 0; j<k-1; j++){
    /// Set values of matrix A, which is if dimension nx(k-1):
    if(j == 0) upper1 += S.getTeams()[j]-1;
    else upper1 = upper2;
    for(int i = lower1; i<= upper1;i++){
      int t = S.getPlayers()[i]-1;
      A(t,j) = (2.0/(S.getTeams()[j] + S.getTeams()[j+1]));
    }
    lower2 = upper1+1;
    upper2 = upper1 + S.getTeams()[j+1];
    for(int i = lower2; i<= upper2;i++){
      int t = S.getPlayers()[i]-1;
      A(t,j) = (2.0/((-1.0)*(S.getTeams()[j] + S.getTeams()[j+1])));
    }
    lower1 = lower2;
    /// Set values of the vectors alpha and beta.
    
  }
  cout << "\nTeamPlayer matrix A:" << endl;
  cout << A << endl;
  
  /// Set vector u and matrix C which are parameters t approximate truncated gaussian with. 
  VectorXd u = A.transpose() * S.getMu();
  VectorXd h = S.getSigma();
  cout << "Sigma Vec: " << endl;
  cout << h << endl;
  for(int i = 0; i<n;i++){
    h[i] +=  S.getBeta() * S.getBeta();
  }
cout << h << endl;
  cout << "\nVector u to approximate truncated Gaussian:" << endl;
  cout << u << endl;
  cout << "\nMatrix C as covariance matrix to approximate truncated Gaussian:" << endl;
  MatrixXd C =  A.transpose() * h.asDiagonal() * A;
  cout << C << endl;
}

/// Test environment for gameSetup. 
int main(void){
  vector<int> p;
  p.push_back(1);
  p.push_back(2);
  p.push_back(3);
  p.push_back(4);
  vector<int> t;
  t.push_back(1);
  t.push_back(2);
  t.push_back(1);
  //t.push_back(1);
  gameSetup S1(p,t);
  S1.printSetup();
  updateAlgortithm(S1,0.1);
}
