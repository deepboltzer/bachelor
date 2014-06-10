#include <math.h>
#include <eigen3/Eigen/Dense>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "gameSetup.h"


using namespace std;
 using namespace Eigen;

/// Players holds the player indices sorted by teams. Teams holds the number of players in team i. 
gameSetup::gameSetup(vector<int> players,vector<int> team)
: sortedPlayers(players),teams(team), numberPlayers(players.size()), numberTeams(team.size()), beta(4.166666666666667), tau(0.08333333333333334){
  this->mu =  VectorXd(this->numberPlayers);
  for(int i = 0; i< numberPlayers;i++){
    this->mu[i] = 25.000; 
  }
  this->sigma = VectorXd(this->numberPlayers);
   for(int i = 0; i< numberPlayers;i++){
    this->sigma[i] = 8.333; 
  }
}

/// Resizes a Team.
void gameSetup::resizeTeam(int teamIndex, int size){
  this->teams.at(teamIndex) = size;
}

void gameSetup::printSetup(){
  cout << "-- This is GameSetup for your created Game: -- \n" << " - Constants are set by recommended default values of Trueskill - \n" << endl;
  cout << "teamVec: \n" << endl;
  
  for(int i = 0; i< this->numberTeams;i++){
    cout << "[" << this->teams.at(i) << "]" << endl;
  }
  cout << " \nplayerVec: \n" << endl;
  
  for(int i = 0; i< this->numberPlayers;i++){
    cout << "[" << this->sortedPlayers.at(i) << "]" << endl;
  }
  
  cout << "\nMu: " << this->mu[0] << endl;
  
  cout << "Sigma: " << this->sigma[0]  << endl;
  
  cout << "Beta: " << this->beta  << endl;
 
  cout << "Tau: " << this->tau << endl;
}

/// Test environment for gameSetup. 
int main(void){
  vector<int> p;
  p.push_back(1);
   p.push_back(2);
    p.push_back(3);
  vector<int> t;
  t.push_back(1);
  t.push_back(2);
  gameSetup S1(p,t);
  S1.printSetup();
}