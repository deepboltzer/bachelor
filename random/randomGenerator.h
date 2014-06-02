#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class randomGenerator{

  
  int sampleIndex;
  int dimension;
  std::vector<double> w;
  
public:
  
  /// Standard Constructor for randomGenerator.
  randomGenerator();
  /// Generate randomGenerator with n samples and sampleIndex.
  randomGenerator(int n, int sampleIndex);
  
  /// Fill Vector w with quasi random Samples.
  void generate_random_sample();
  
};
