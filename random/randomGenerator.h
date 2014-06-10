#include <iostream>
#include <cmath>
#include <vector>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

class randomGenerator{

  
  int sampleIndex;
  int dimension;
  VectorXd w;
  
public:
  
  /// Standard Constructor for randomGenerator.
  randomGenerator();
  /// Generate randomGenerator with n samples and sampleIndex.
  randomGenerator(int n, int sampleIndex);
  
  /// Fill Vector w with quasi random Samples.
  void generateRandomSample();
  
  /// PrintRandomSample.
  void printSample();
  
};
