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

/// Returns number of players of the setup.
double gameSetup::getNumberPlayers(){
  return this->numberPlayers;
}

/// Returns number of teams of the setup.
double gameSetup::getNumberTeams(){
  return this->numberTeams;
}

/// Returns tau;
double gameSetup::getTau(){
  return this-> tau;
}

/// Returns beta;
double gameSetup::getBeta(){
  return this-> beta;
}

/// Returns variance vector.
VectorXd gameSetup::getSigma(){
  return this->sigma;
}

/// Returns mean mu as a vector. 
VectorXd gameSetup::getMu(){
  return this->mu;
}

/// Returns team vector.
vector<int> gameSetup::getTeams(){
  return this->teams;
}

/// Returns sortedPlayer vector.
vector<int> gameSetup::getPlayers(){
  return this->sortedPlayers;
}
/// Returns variance for a player with playerIndex.
double gameSetup::getSigmaByIndex(int playerIndex){
  return this->sigma[playerIndex];
}
/// Set the varaince player with playerIndex to value.
void gameSetup::setSigma(int playerIndex, double value){
  this->sigma[playerIndex] = value;
}

// Prints gameSetup.
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

