#include <math.h>
#include <eigen3/Eigen/Dense>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "gameSetup.h"

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
  MatrixXd A(n,k);
  
  /// Set values of matrix A:
  bool first = true;
  int lower1 = 0;
  int upper1 = 0;
  int lower2 = 0;
  int upper2 = 0;
  
  for(int j = 0; j<k-1; j++){
    /// Set values of matrix A:
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
  for(int i = 0; i<n;i++){
    h[i] +=  S.getBeta() * S.getBeta();
  }
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
  p.push_back(5);
  vector<int> t;
  t.push_back(2);
  t.push_back(2);
  t.push_back(1);
  //t.push_back(1);
  gameSetup S1(p,t);
  S1.printSetup();
  updateAlgortithm(S1,0.1);
}