#include <iostream>

/// A Gaussian distribution based on float numbers.
/// in exponential parameterisation. 
class Gaussian{
  float PrecisionMean;
  float Precision;
  float Mu;
  float Mean; 
  float Variance;
  float StadardDeviation;
  float Sigma;
};

/// Creates a Gaussian from Precision and PrecisionMean. 
Gaussian::Gaussian (float PrecisionMean, float Precision) {
  
  this.PrecisionMean = PrecisionMean;
  this.Precision = Precision;
  this.Mu = PrecisionMean / Precision;
  this.Mean = this.Mu;
  this.Variance = 1./Precision;
  this.StadardDeviation = sqrt(this.Variance);
  this.Sigma = this.StandardDeviation;
}

/// Creates a Gaussian in (mean,standard-deviation) coordinates
Gaussian::Gaussian (float mu, float sigma){
  float sigmaSquared = sigma*sigma;
  new Gaussian(mu/(sigmaSquared),1./(sigmaSquared));
}

/// Multiplies two Gaussians  
static Gaussian multiply_Gaussian(Gaussian a, Gaussian b){
  return new Gaussian(a.PrecisionMean + b.PrecisionMean, a.Precision + b.Precision);
}

/// Divides two Gaussians
static Gaussian divide_Gaussian(Gaussian a, Gaussian b){
  return new Gaussian(a.PrecisionMean - b.PrecisionMean, a.Precision - b.Precision);
}

/// Computes the absolute difference between two Gaussians
static float AbsoluteDifference(Gaussian a, Gaussian b){
  return max(abs(a.PrecisionMean - b.PrecisionMean),sqrt(abs(a.Precision - b.Precision)));
}        

/// Computes the log-normalisation factor when two normalised Gaussians gets multiplied
static float LogProductNormalisation(Gaussian a, Gaussian b){
  if (a.Precision == 0.0) return 0.0;
  else if (b.Precision == 0.0) return 0.0; 
  else {
    float varSum = a.Variance + b.Variance;
    float muDiff = a.Mean - b.Mean;
    return -0.91893853320467267 - log(varSum)/2.0 - muDiff*muDiff/(2.0 * varSum);
  }
}
                
 /// Computes the log-normalisation factor when two normalised Gaussians gets divided
static float LogRatioNormalisation (Gaussian a, Gaussian b){
  if (a.Precision == 0.0) return 0.0;
  else if (b.Precision == 0.0) return 0.0;
  else{
    float v2 = b.Variance;
    float varDiff = v2 - a.Variance;
    float muDiff = a.Mean - b.Mean;
    if(varDiff == 0.0) return 0.0;
    else return log(v2) + 0.91893853320467267 - log(varDiff)/2.0 + muDiff*muDiff/(2.0 * varDiff);
  }
}

/// Computes the complementary error function. This function is defined     
/// by 2/sqrt(pi) * integral from x to infinity of exp (-t^2) dt.
static float erfc(float x){
  if(IsNegativeInfinity(x) == true) return 2.0;
  else if(IsPositiveInfinity(x) == true) return 0.0;
  else{
    float z = abs(x);
    float t = 1.0/(1.0 + 0.5 * z);
    float res = t * exp (-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    if (x >= 0.0) return res;
    else return 2.0 - res;
  }
}

/// Computes the inverse of the complementary error function
    
static float erfcinv(float y){
  if (y < 0.0 || y > 2.0) return 0;
  else if(y == 0.0) return PositiveInfinity;
  else if(y == 2.0) return NegativeInfinity;
  else{
    float x;
    if (y >= 0.0485 && y <= 1.9515){
                    float q = y - 1.0;
                    float r = q * q;
                     x = (((((0.01370600482778535*r - 0.3051415712357203)*r + 1.524304069216834)*r - 3.057303267970988)*r + 2.710410832036097)*r - 0.8862269264526915) * q /
                    (((((-0.05319931523264068*r + 0.6311946752267222)*r - 2.432796560310728)*r + 4.175081992982483)*r - 3.320170388221430)*r + 1.0);
                else if (y < 0.0485) 
                    float q = sqrt (-2.0 * log (y / 2.0));
                    x = (((((0.005504751339936943*q + 0.2279687217114118)*q + 1.697592457770869)*q + 1.802933168781950)*q + -3.093354679843504)*q - 2.077595676404383) /
                    ((((0.007784695709041462*q + 0.3224671290700398)*q + 2.445134137142996)*q + 3.754408661907416)*q + 1.0);
                else if (y > 1.9515) 
                    float q = sqrt (-2.0 * log (1.0 - y / 2.0));
                    x = (-(((((0.005504751339936943*q + 0.2279687217114118)*q + 1.697592457770869)*q + 1.802933168781950)*q + -3.093354679843504)*q - 2.077595676404383) /
                     ((((0.007784695709041462*q + 0.3224671290700398)*q + 2.445134137142996)*q + 3.754408661907416)*q + 1.0));
                else x = 0.0;
    }
            float u = (erfc (x) - y) / (-2.0 / sqrt (Math.PI) * exp (-x * x));
            return x - u / (1.0 + x * u);
  }
}
/// Phi   
/// Computes the cummulative Gaussian distribution at a specified point of interest
static float normcdf(float t) {
  float sqrt2 = 1.4142135623730951; 
  return ((erfc (-t / sqrt2)) / 2.0);
}

/// Computes the Gaussian density at a specified point of interest
static float normpdf (float t){
  float invsqrt2pi = 0.398942280401433;
  return (invsqrt2pi * exp (- (t * t / 2.0)));
} 

/// PhiInverse
/// Computes the inverse of the cummulative Gaussian distribution (qunatile function) at a specified point of interest
static float norminv(float p){ 
  float sqrt2 = 1.4142135623730951; 
  return (-sqrt2 * erfcinv (2.0 * p));
}


/// Computes the additive correction of a single-sided truncated Gaussian with unit variance
static float v_t_epsilon(float t, float epsilon){
  float valnormcdf = normcdf(t-epsilon);
  if(valnormcdf < 2.222758749e-162) return (-t+epsilon);
  else return (normpdf(t - epsilon) / valnormcdf); 
} 

/// Computes the multiplicative correction of a single-sided truncated Gaussian with unit variance
static float w_t_epsilon(float t, float epsilon){
  float valnormcdf = normcdf(t-epsilon);
  if(valnormcdf < 2.222758749e-162){
    if (t < 0.0) return 1.0;
    else return 0.0;
  }
  else{
    float vt = v_t_epsilon(t, epsilon);
    return vt * (vt + t - epsilon);
  }
}

/// Computes the additive correction of a double-sided truncated Gaussian with unit variance
static float v0_t_epsilon(float t, float epsilon){
  float valnormcdf = normcdf(t-epsilon) - normcdf(-epsilon - v);
  if(valnormcdf  < 2.222758749e-162){
    if(t < 0.0) return (-t-epsilon);
    else return ((-t)+epsilon);
  }
  else {
   float num = normpdf(-epsilon-v) - normpdf(epsilon -v);
   if(t < 0.0) return (-num/valnormcdf); 
   else return (num/valnormcdf);
  } 
}  

/// Computes the multiplicative correction of a double-sided truncated Gaussian with unit variance
static float w0_t_epsilon(float t, float epsilon){
  float v = abs(t);
  valnormcdf = normcdf(epsilon - v) - normcdf(-epsilon - v);
  if(valnormcdf  < 2.222758749e-162) return 1.0;
  else{
    float vt = v0_t_epsilon(v,epsilon);
    return vt*vt + ((epsilon-v) * normpdf (epsilon-v) - (-epsilon-v) * normpdf (-epsilon-v))/valnormcdf;
  }
}

int main()
{
  std::cout << "Hello World!";
}