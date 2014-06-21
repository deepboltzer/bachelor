#include <math.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Cholesky>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "gameSetup.h"
#include <boost/math/distributions/normal.hpp>
#include "../gaussians/gaussian.h"
#include "../random/randomGenerator.h"

using namespace std;
double ninf  = -numeric_limits<double>::infinity();
double inf  = numeric_limits<double>::infinity();
void genzAlgorithm(VectorXd a, VectorXd b, MatrixXd covar) 
{
  boost::math::normal norm;
  MatrixXd C(3,3); 
  C << 1.0 , 0., 0. ,
  3./5., 4./5., 0.,
  1./3., 2./3.,2./3.;
  int N, m;
  double varsum, itsum;
  N = 0; varsum = 0; itsum = 0; m = 3; 
  double d,e;
  if(a[0] == ninf) d = 0;
  else d = cdf(norm,(a[1]/C(1,1)));
  if(b[0] = inf) e = 1; 
  else e = cdf(norm,(b[0]/C(1,1)));
  double f = e-d; 
  cout << "f: " << f << endl;
  double delta = 0;
  double error = 1000;
  VectorXd y(m);
  double sum = 0; 
  cout << d << "	"<< e << endl;
  randomGenerator rg(m,1);
  rg.generateRandomSample();
  VectorXd w = rg.getSample();
  cout << w << endl;
  while(N<=1000){
  for(int i = 1; i<m; i++){
    w[0] = 0.47839;
    sum = 0; 
    cout << "w: "<< w[i-1] << "d + w[i-1]*(e-d): " << d + w[i-1]*(e-d) << endl;
    y[i-1] = quantile(norm,(d + w[i-1]*(e-d))); 
    cout << "y: " << y[i-1] << endl;
    for(int j = 0; j<= i-1;j++){
      sum += C(i,j)*y[j];
    }
    cout << "sum:" << sum << endl;
    if(a[i] == ninf) d = 0;
    else d = cdf(norm,((a[i] -sum)/C(i,i))); 
    if(b[i] == inf) e = 1; 
    else e = cdf(norm,((b[i] -sum)/C(i,i)));
    cout << d << "	"<< e << endl;
    f = (e - d)*f;
    cout << "(b[i] -sum)/C(i,i)): " << (b[i] -sum)<< endl;
  }
  N = N+1;
  delta = (f - itsum)/double(N);
  itsum = itsum + delta;
  varsum = double(N-2)*varsum/double(N) + delta*delta;
  error = 2.5*sqrt(varsum);
  cout << "intsum: " << itsum << endl;
  cout << "error: " << error << endl;
  }
}
int main(void){
  VectorXd a(3);a << ninf,ninf,ninf;
  VectorXd b(3);b << 1,4,2;
  MatrixXd covar(3,3); covar <<  1.0, 3./5.,       1./3.,
				 5./3.,  1.,      11./15.,
				  1./3.      ,  11./15.      ,       1.; 
  genzAlgorithm(a,b,covar);
  boost::math::normal norm;
}