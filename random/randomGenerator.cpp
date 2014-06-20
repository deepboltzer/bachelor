#include <iostream>
#include "randomGenerator.h"
#include <cmath>
#include <vector>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

/// Generate randomGenerator with n samples and sampleIndex.
randomGenerator::randomGenerator(int n, int j)
: dimension(n), sampleIndex(j)
{ 
  this->w = VectorXd(this->dimension);
}

/// Fill Vector w with quasi random Samples.
void randomGenerator::generateRandomSample()
{
    int i; 
    for(i=0;i<this->dimension;i++){
      this->w(i) = std::abs(2*((this->sampleIndex)*pow(2,double(i)/double(this->dimension+1))- floor((this->sampleIndex)*pow(2,double(i)/double(this->dimension +1)))) - 1);
    }
}

VectorXd randomGenerator::getSample()
{
    return this->w;
}
/// Print RandomSample.
void randomGenerator::printSample(){
  for(int i=0;i<this->dimension;i++){
    std::cout << this->w(i) << std::endl;
  }
}


