#include <math.h>
#include <eigen3/Eigen/Dense>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "gameSetup.h"

/// Generate gameSetup from number of Players.
using namespace std;

gameSetup::gameSetup(int p, int t)
: numberPlayers(p), numberTeams(t), beta(4.166666666666667), tau(0.08333333333333334) 
{
  /// Set mu and sigma to default trueskill values. 
  this->mu =  vector<double>(this->numberPlayers,25.000);
  this->sigma = vector<double>(this->numberPlayers,8.333);
  this->teams = vector<int> (t);
}

/// Add player to team. 
void gameSetup::addPlayerToTeam(int teamIndex){
  this->numberPlayers++;
  this->teams[teamIndex-1]++;
  vector<int>::iterator iter;
  int sum = 0;
  for(iter = teams.begin(); iter != teams.end(); ++iter){
    sum = sum + (*iter);
  }
  iter = sortedPlayers.begin() + sum;
  this->sortedPlayers.insert(iter,this->numberPlayers);
}

void gameSetup::printSetup(){
  cout << "teamVec: \n" << endl;
  
  for(int i = 0; i< this->numberTeams;i++){
    cout << "[" << this->teams.at(i) << "]" << endl;
  }
  cout << " \n playerVec: \n" << endl;
  
  for(int i = 0; i< this->numberPlayers;i++){
    cout << "[" << this->sortedPlayers.at(i) << "]" << endl;
  }
  
  cout << "\nMu: " << this->mu.at(0) << "\n" << endl;
  
  cout << "Sigma: " << this->sigma.at(0) << "\n"  << endl;
  
  cout << "Beta: " << this->beta << "\n"  << endl;
 
  cout << "Tau: " << this->tau << "\n" << endl;
}

int main(void){
  
}