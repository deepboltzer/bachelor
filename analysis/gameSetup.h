#include <math.h>
#include <eigen3/Eigen/Dense>
#include <stdlib.h>
#include <iostream>
#include <vector>


using namespace std;
using namespace Eigen;

class gameSetup{
  
public:
  int numberPlayers;
  int numberTeams;
  /// Holds player indices sorted by teams.
  vector<int> sortedPlayers;
  /// Holds numbers of players sorted by teams.
  vector<int> teams;
  VectorXd mu;
  VectorXd sigma;
  double beta;
  double tau;

  /// Default Constructor for gameSetup;
  gameSetup();
  
  /// Generate gameSetup from number of Players.
  gameSetup(vector<int> players,vector<int> team);
  
  /// Returns number of players of the setup.
  double getNumberPlayers();

  /// Returns number of players of the setup.
  double getNumberTeams();

  /// Returns tau;
  double getTau();

  /// Returns Variance for a player with playerIndex. 
  double getSigma(int playerIndex);

  /// Set the varaince player with playerIndex to value.
  void setSigma(int playerIndex, double value);
  
  /// Returns team vector.
  vector<int> getTeams();
  
  /// Returns team vector.
  vector<int> getPlayers();
  
  /// Resizes a team.
  void resizeTeam(int teamIndex, int size);

  
  /// Prints Setup. 
  void printSetup();
};