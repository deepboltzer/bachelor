class Gaussian{
 public:
  double PrecisionMean;
  double Precision;
  double Mu;
  double Mean; 
  double Variance;
  double StadardDeviation;
  double Sigma;
  

  /// Default Constructor.
  Gaussian();
  /// Constructor which generates Gaussian from Precision and PrecisionMean.
  Gaussian(double PrecisionMean, double Precision);
  /// Generate Gaussian from mu ans sigma. 
  Gaussian(double mu, double sigma);
  /// Destructor for Gaussians.
  ~Gaussian();
};