#include <iostream>
#include "gaussian.h"
#include <cmath>
#include <limits>
#include <fstream> 
#include "../gnuplot-iostream/gnuplot-iostream.h"

using namespace std;

// Demo of vector plot.
// Compile it with:
//   g++ -o example-vector example-vector.cc -lboost_iostreams -lboost_system -lboost_filesystem
/// Test environment for Gaussians:
int main(void){
  Gnuplot gp;
  Gaussian G1(0,1);
  int i;
  std::fstream f1,f2;
  f1.open("cdf_01gaussian.dat", ios::out);
  f2.open("pdf_01gaussian.dat", ios::out);
  //for(i=-500;i<=500;i++){
    // f1 << G1.normcdf(double(i)/100.0) << "\n" << endl;
     //f2 << G1.normpdf(double(i)/100.0) << "\n" << endl;
  //}
  f1.close();
  f2.close();
  //Gaussian G2(0.25,0.75);
  //Gaussian* G3 = (G1.multiply_Gaussian(G2));
  //cout << G1.normpdf(0) << endl;
}