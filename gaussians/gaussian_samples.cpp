#include <iostream>
#include "gaussian.h"
#include "truncated_gaussian.h"
#include <cmath>
#include <limits>
#include <fstream> 
#include "../gnuplot-iostream/gnuplot-iostream.h"

using namespace std;


/// Test environment for Gaussians:
/// Compile it with:
/// g++ -Wall -o gaussian gaussian.cpp gaussian_samples.cpp -lboost_iostreams -lboost_system -lboost_filesystem
/*
int plotGaussianDensity(){
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
		double y = G1.at(double(i)/100, 0,1);
		y_pts.push_back(y);
		double x = G2.at(double(i)/100,0,2);
		x_pts.push_back(x);
		double z = G3.at(double(i)/100, -1, 2);
		z_pts.push_back(z);
	}

	std::cout << "Creating Gaussians.png" << std::endl;
	gp << "set output 'Gaussians.png'\n";
	//gp << "set xrange [-5:5]\n";
	gp << "set xlabel 't'\n";
	gp << "plot '-' with lines title 'N(0,1)', '-' with lines title 'N(0,4)', '-' with lines title 'N(-1,4)'\n";
	gp.send1d(y_pts);
	gp.send1d(x_pts);
	gp.send1d(z_pts);
  return 0; 
}*/
int main(void){
 Gnuplot gp;
truncatedGaussianCorrectionFunctions tg;
	std::vector<double> mu_corr1;
	std::vector<double> mu_corr2;
	std::vector<double> mu_corr3;
	std::vector<double> corr_func;
	for(int i=-500;i<10;i++) {
	  double corr1 = tg.vWithinMargin(double(i)/100, 0.1);
	  mu_corr1.push_back(corr1);
	}
	for(int i=-500;i<200;i++) {
	  double corr2 = tg.vWithinMargin(double(i)/100, 2.0);
	  mu_corr2.push_back(corr2);
	} 
	for(int i=-500;i<400;i++) {
	  double corr3 = tg.vWithinMargin(double(i)/100, 4.0);
	  mu_corr3.push_back(corr3);
	} 
	
	
	for(int i=10;i<=500;i++) {
	  double corr1 = tg.vExceedsMargin(double(i)/100, 0.1);
	  mu_corr1.push_back(corr1);
	}
	for(int i=200;i<=500;i++) {
	  double corr2 = tg.vExceedsMargin(double(i)/100, 2.0);
	  mu_corr2.push_back(corr2);
	} 
	for(int i=400;i<=500;i++) {
	  double corr3 = tg.vExceedsMargin(double(i)/100, 4.0);
	  mu_corr3.push_back(corr3);
	} 
	for(int i=-500;i<=500;i++) {
	  corr_func.push_back(tg.at(double(i)/100));
	}
	
	std::cout << "Creating SigmaCorrection.png" << std::endl;
	gp << "set output 'SigmaCorrection.png'\n";
	gp << "set grid\n";
	gp << "plot '-' with points title 'epsilon = 0.50', '-' with points title 'epsilon = 1.00', '-' with points title 'epsilon = 4.00', '-' with points title 'truncated error func'\n";
	//gp << "set grid \n";
	//gp << "set xlabel 't', set ylabel 'v(t,a,b)'\n";
	gp.send1d(mu_corr1);
	gp.send1d(mu_corr2);
	gp.send1d(mu_corr3);
	gp.send1d(corr_func);

}
