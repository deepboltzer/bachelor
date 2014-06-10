#include <math.h>
#include <eigen3/Eigen/Dense>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "gameSetup.h"

void updateAlgortithm(gameSetup S, double drawMargin){
  int n = S.getNumberPlayers();
  int k = S.getNumberTeams();
  /// First Update of player variance:
  for(int j = 0; j< k; j++){
    double value = sqrt(S.getSigma(j)*S.getSigma(j) + S.getTau()*S.getTau());
    S.setSigma(j,value);
  }
  
  /// Declaration of PlayerTeamMatrix A; 
  MatrixXd A(n,k);
  
  ///
  bool first = true;
  int lower1 = 0;
  int upper1 = 0;
  int lower2 = 0;
  int upper2 = 0;
  
  for(int j = 0; j<k-1; j++){
    cout << "j:" << j << endl;
    if(j == 0) upper1 += S.getTeams()[j]-1;
    else upper1 = upper2;
     cout << "upper1:" << upper1 << endl;
    for(int i = lower1; i<= upper1;i++){
      int t = S.getPlayers()[i]-1;
      cout << "t1:" << S.getTeams()[j] << endl;
      A(t,j) = (2.0/(S.getTeams()[j] + S.getTeams()[j+1]));
    }
    lower2 = upper1+1;
    cout << "lower1:" << lower2 << endl;
    upper2 = upper1 + S.getTeams()[j+1];
    cout << "upper2:" << upper2 << endl;
    for(int i = lower2; i<= upper2;i++){
      int t = S.getPlayers()[i]-1;
      cout << "t2:" << t << endl;
      A(t,j) = (2.0/((-1.0)*(S.getTeams()[j] + S.getTeams()[j+1])));
    }
    lower1 = lower2;
    cout << "lower2:" << lower1 << endl;
    
  }
  std::cout << A << endl;
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
  t.push_back(1);
  t.push_back(2);
  t.push_back(1);
  t.push_back(1);
  gameSetup S1(p,t);
  S1.printSetup();
  updateAlgortithm(S1,0.1);
}