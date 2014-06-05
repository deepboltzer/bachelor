#include <iostream>
#include "gaussian.h"
#include <cmath>
#include <limits>
#include <fstream> 
#include "../gnuplot-iostream/gnuplot-iostream.h"

using namespace std;


/// Test environment for Gaussians:
/// Compile it with:
/// g++ -Wall -o gaussian gaussian.cpp gaussian_samples.cpp -lboost_iostreams -lboost_system -lboost_filesystem

int main(void){
 Gnuplot gp;
  Gaussian G1(0,1);
  Gaussian G2(0,4);
  Gaussian G3(-1,4);
  //cout << G1.normpdf(0) << endl;
  gp << "set terminal png\n";

	std::vector<double> y_pts;
	std::vector<double> x_pts;
	std::vector<double> z_pts;
	for(int i=-500;i<=500;i++) {
		double y = G1.normpdf(double(i)/100.0);
		y_pts.push_back(y);
		double x = G2.normpdf(double(i)/100.0);
		x_pts.push_back(x);
		double z = G3.normpdf(double(i)/100.0);
		z_pts.push_back(z);
	}

	std::cout << "Creating Gaussians.png" << std::endl;
	gp << "set output 'Gaussians.png'\n";
	gp << "plot '-' with lines title 'N(0,1)', '-' with lines title 'N(0,4)', '-' with lines title 'N(-1,4)'\n";
	gp.send1d(y_pts);
	gp.send1d(x_pts);
	gp.send1d(z_pts);
}
