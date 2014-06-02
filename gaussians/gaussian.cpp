#include <iostream>
#include "gaussian.h"
#include <cmath>
#include <limits>

using namespace std;

/// Creates a Gaussian in (mean,standard-deviation) coordinates
Gaussian::Gaussian(double mu, double sigma)
: PrecisionMean(mu/(sigma*sigma)), Precision(1./(sigma*sigma)), Mu(mu), Mean(mu), Variance(sigma*sigma), StadardDeviation(sigma), Sigma(sigma) 
{
}

/// Deconstructor for Gaussians

Gaussian::~Gaussian()
{
  
}
/// Multiplies two Gaussians  
Gaussian* Gaussian::multiply_Gaussian(Gaussian a, Gaussian b){
  return new Gaussian(a.PrecisionMean + b.PrecisionMean, a.Precision + b.Precision);
}

/// Divides two Gaussians
Gaussian* Gaussian::divide_Gaussian(Gaussian a, Gaussian b){
  return new Gaussian(a.PrecisionMean - b.PrecisionMean, a.Precision - b.Precision);
}

/// Computes the absolute difference between two Gaussians
double Gaussian::AbsoluteDifference(Gaussian a, Gaussian b){
  return max<double>(std::abs(a.PrecisionMean - b.PrecisionMean),sqrt(std::abs(a.Precision - b.Precision)));
}        

/// Computes the log-normalisation factor when two normalised Gaussians gets multiplied
double Gaussian::LogProductNormalisation(Gaussian a, Gaussian b){
  if (a.Precision == 0.0) return 0.0;
  else if (b.Precision == 0.0) return 0.0; 
  else {
    double varSum = a.Variance + b.Variance;
    double muDiff = a.Mean - b.Mean;
    return (-0.91893853320467267 - log(varSum)/2.0 - muDiff*muDiff/(2.0 * varSum));
  }
}
                
/// Computes the log-normalisation factor when two normalised Gaussians gets divided
double Gaussian::LogRatioNormalisation (Gaussian a, Gaussian b){
  if (a.Precision == 0.0) return 0.0;
  else if (b.Precision == 0.0) return 0.0;
  else{
    double v2 = b.Variance;
    double varDiff = v2 - a.Variance;
    double muDiff = a.Mean - b.Mean;
    if(varDiff == 0.0) return 0.0;
    else return (log(v2) + 0.91893853320467267 - log(varDiff)/2.0 + muDiff*muDiff/(2.0 * varDiff));
  }
}

///Tests if a double is PositiveInfinity

bool Gaussian::IsPositiveInfinity(double x){
  
 if(x == std::numeric_limits<double>::infinity()) return true;
 else return false;
}

/// Tests if a double is NegativeInfinity

bool Gaussian::IsNegativeInfinity(double x){
  if(x == -(std::numeric_limits<double>::infinity())) return true;
  else return false;
}

/// Computes the complementary error function. This function is defined     
/// by 2/sqrt(pi) * integral from x to infinity of exp (-t^2) dt.
double Gaussian::erfc(double x){
  if(IsNegativeInfinity(x) == true) return 2.0;
  else if(IsPositiveInfinity(x) == true) return 0.0;
  else{
    double z = std::abs(x);
    double t = 1.0/(1.0 + 0.5 * z);
    double res = t * exp (-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    if (x >= 0.0) return res;
    else return 2.0 - res;
  }
}

/// Computes the inverse of the complementary error function
    
double Gaussian::erfcinv(double y){
  if (y < 0.0 || y > 2.0) return 0;
  else if(y == 0.0) return std::numeric_limits<double>::infinity();
  else if(y == 2.0) return -(std::numeric_limits<double>::infinity());
  else{
    double x;
    if (y >= 0.0485 && y <= 1.9515){
                    double q = y - 1.0;
                    double r = q * q;
                     x = (((((0.01370600482778535*r - 0.3051415712357203)*r + 1.524304069216834)*r - 3.057303267970988)*r + 2.710410832036097)*r - 0.8862269264526915) * q /
                    (((((-0.05319931523264068*r + 0.6311946752267222)*r - 2.432796560310728)*r + 4.175081992982483)*r - 3.320170388221430)*r + 1.0);
    }
    else if (y < 0.0485){ 
                    double q = sqrt (-2.0 * log (y / 2.0));
                    x = (((((0.005504751339936943*q + 0.2279687217114118)*q + 1.697592457770869)*q + 1.802933168781950)*q + -3.093354679843504)*q - 2.077595676404383) /
                    ((((0.007784695709041462*q + 0.3224671290700398)*q + 2.445134137142996)*q + 3.754408661907416)*q + 1.0);
    }
    else if (y > 1.9515){
                    double q = sqrt (-2.0 * log (1.0 - y / 2.0));
                    x = (-(((((0.005504751339936943*q + 0.2279687217114118)*q + 1.697592457770869)*q + 1.802933168781950)*q + -3.093354679843504)*q - 2.077595676404383) /
                     ((((0.007784695709041462*q + 0.3224671290700398)*q + 2.445134137142996)*q + 3.754408661907416)*q + 1.0));
    }
    else x = 0.0;
    double u = (erfc (x) - y) / (-2.0 / sqrt (M_PI) * exp (-x * x));
    return x - u / (1.0 + x * u);
  }
}
/// Phi   
/// Computes the cummulative Gaussian distribution at a specified point of interest
double Gaussian::normcdf(double t) {
  double sqrt2 = 1.4142135623730951; 
  return ((erfc (-t / sqrt2)) / 2.0);
}

/// Computes the Gaussian density at a specified point of interest
double Gaussian::normpdf (double t){
  double invsqrt2pi = 0.398942280401433;
  return (invsqrt2pi * exp (- (t * t / 2.0)));
} 

/// PhiInverse
/// Computes the inverse of the cummulative Gaussian distribution (qunatile function) at a specified point of interest
double Gaussian::norminv(double p){ 
  double sqrt2 = 1.4142135623730951; 
  return (-sqrt2 * erfcinv (2.0 * p));
}


/// Computes the additive correction of a single-sided truncated Gaussian with unit variance
double Gaussian::v_t_epsilon(double t, double epsilon){
  double valnormcdf = normcdf(t-epsilon);
  if(valnormcdf < 2.222758749e-162) return (-t+epsilon);
  else return (normpdf(t - epsilon) / valnormcdf); 
} 

/// Computes the multiplicative correction of a single-sided truncated Gaussian with unit variance
double Gaussian::w_t_epsilon(double t, double epsilon){
  double valnormcdf = normcdf(t-epsilon);
  if(valnormcdf < 2.222758749e-162){
    if (t < 0.0) return 1.0;
    else return 0.0;
  }
  else{
    double vt = v_t_epsilon(t, epsilon);
    return vt * (vt + t - epsilon);
  }
}

/// Computes the additive correction of a double-sided truncated Gaussian with unit variance
double Gaussian::v0_t_epsilon(double t, double epsilon){
  double v = std::abs(t);
  double valnormcdf = normcdf(t-epsilon) - normcdf(-epsilon - v);
  if(valnormcdf  < 2.222758749e-162){
    if(t < 0.0) return (-t-epsilon);
    else return ((-t)+epsilon);
  }
  else {
   double num = normpdf(-epsilon-v) - normpdf(epsilon -v);
   if(t < 0.0) return (-num/valnormcdf); 
   else return (num/valnormcdf);
  } 
}  

/// Computes the multiplicative correction of a double-sided truncated Gaussian with unit variance
double Gaussian::w0_t_epsilon(double t, double epsilon){
  double v = std::abs(t);
  double valnormcdf = normcdf(epsilon - v) - normcdf(-epsilon - v);
  if(valnormcdf  < 2.222758749e-162) return 1.0;
  else{
    double vt = v0_t_epsilon(v,epsilon);
    return vt*vt + ((epsilon-v) * normpdf (epsilon-v) - (-epsilon-v) * normpdf (-epsilon-v))/valnormcdf;
  }
}

int main(void){
  
}