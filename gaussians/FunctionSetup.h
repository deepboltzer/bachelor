
using namespace Eigen;

class FunctionSetup{
private: 
 VectorXd mean;
 MatrixXd covar;
 MatrixXd A;
 double beta;

public: 
  /// Default Constructor.
  FunctionSetup();
  /// Constructor which generates FunctionSetup with param mean,CovarianceMatrix, Matrix A and trueskill factor beta.
  FunctionSetup(VectorXd m, MatrixXd S, MatrixXd partial, double b);
  /// Generates u(mu) = A^T*mu;
  VectorXd function_u();
};