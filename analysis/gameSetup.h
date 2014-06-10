#include <math.h>
#include <eigen3/Eigen/Dense>
#include <stdlib.h>
#include <iostream>
#include <vector>


using namespace std;

class gameSetup{
  int numberPlayers;
  int numberTeams;
  /// Holds player indices sorted by teams.
  vector<int> sortedPlayers;
  /// Holds numbers of players sorted by teams.
  vector<int> teams;
  vector<double> mu;
  vector<double> sigma;
  double beta;
  double tau;

public:
  /// Default Constructor for gameSetup;
  gameSetup();
  /// Generate gameSetup from number of Players.
  gameSetup(int p, int t);
  
  /// Add player to team. 
  void addPlayerToTeam(int teamIndex);
  
  /// Prints Setup. 
  void printSetup();
};