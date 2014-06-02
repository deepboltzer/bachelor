#include <iostream>
#include "randomGenerator.h"
#include <cmath>
#include <vector>

using namespace std;

/// Generate randomGenerator with n samples and sampleIndex.
randomGenerator::randomGenerator(int n, int j)
: dimension(n), sampleIndex(j)
{ 
  std::vector<double> w (n);
}

/// Fill Vector w with quasi random Samples.
void randomGenerator::generate_random_sample()
{
    int i; 
    for(i=0;i<this->dimension;i++){
      this->w.at(i) = std::abs(2*((this->sampleIndex)*pow(2,double(i)/double(this->dimension+1))- floor((this->sampleIndex)*pow(2,double(i)/double(this->dimension +1)))) - 1);
    }
}

int main(void){
  
}