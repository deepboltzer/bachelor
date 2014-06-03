
/// A Gaussian distribution based on double numbers.
/// in exponential parameterisation. 

class Gaussian{
private:
  double PrecisionMean;
  double Precision;
  double Mu;
  double Mean; 
  double Variance;
  double StadardDeviation;
  double Sigma;
  
 public: 

  /// Default Constructor.
  Gaussian();
  /// Generate Gaussian from mu ans sigma. 
  Gaussian(double mu, double sigma);
  /// Destructor for Gaussians.
  ~Gaussian();
  
  /// Multiplies two Gaussians  
  Gaussian* multiply_Gaussian(Gaussian b);
  /// Divides two Gaussians
  Gaussian* divide_Gaussian(Gaussian b);
  /// Computes the absolute difference between two Gaussians
  double AbsoluteDifference(Gaussian a, Gaussian b);
  /// Computes the log-normalisation factor when two normalised Gaussians gets multiplied
  double LogProductNormalisation(Gaussian a, Gaussian b);
  /// Computes the log-normalisation factor when two normalised Gaussians gets divided
  double LogRatioNormalisation (Gaussian a, Gaussian b);
  ///Tests if a double is PositiveInfinity
  bool IsPositiveInfinity(double x);
  /// Tests if a double is NegativeInfinity
  bool IsNegativeInfinity(double x);
  /// Computes the complementary error function.
  double erfc(double x);
  /// Computes the inverse of the complementary error function
  double erfcinv(double y);
  /// Computes the cummulative Gaussian distribution at a specified point of interest
  double normcdf(double t);
  /// Computes the Standard Gaussian density at a specified point of interest
  double standard_normpdf (double t);
  /// Computes the Gaussian density at a specified point of interest
  double normpdf(double t);
  /// Computes the inverse of the cummulative Gaussian distribution (qunatile function) at a specified point of interest
  double norminv(double p);
  /// Computes the additive correction of a single-sided truncated Gaussian with unit variance
  double v_t_epsilon(double t, double epsilon);
  /// Computes the multiplicative correction of a single-sided truncated Gaussian with unit variance
  double w_t_epsilon(double t, double epsilon);
  /// Computes the additive correction of a double-sided truncated Gaussian with unit variance
  double v0_t_epsilon(double t, double epsilon);
  /// Computes the multiplicative correction of a double-sided truncated Gaussian with unit variance
  double w0_t_epsilon(double t, double epsilon);
  
};