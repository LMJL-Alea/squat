#include "rotations.h"

Eigen::Matrix3d logSO3C(const Eigen::MatrixXd &R)
{
  double workValue = 0.5 * R.trace() - 0.5;
  if (workValue > 1.0)
    workValue = 1.0;
  if (workValue < -1.0)
    workValue = -1.0;
  double theta = std::acos(workValue);
  double denomValue = 2.0 * std::sin(theta);

  if (denomValue < std::sqrt(std::numeric_limits<double>::epsilon()))
    return Eigen::Matrix3d::Zero();

  return (R - R.transpose()) * theta / denomValue;
}

Eigen::Matrix3d expskewC(const Eigen::MatrixXd &M)
{
  /* This function takes a 3-by-3 skew symmetric matrix (in so(3)) and
   returns the exponential, a 3-by-3 rotations (in SO(3)) */

  double MMt = (M - M.transpose()).sum();

  if (std::abs(MMt) > 0.01)
    throw Rcpp::exception("The expskewC function is expecting a 3-by-3 skew symmetric matrix.");

  Eigen::Matrix3d expM = Eigen::Matrix3d::Identity();

  double a = std::sqrt(0.5 * (M.transpose() * M).trace());

  if (std::abs(a) < 1.0e-6)
    return expM;

  expM += (std::sin(a) / a) * M + (1.0 - std::cos(a)) * std::pow(a, -2.0) * M * M;

  return expM;
}

Eigen::Matrix3d projectSO3C(const Eigen::MatrixXd &M)
{
  /* This function will project an arbitrary 3-by-3 matrix M in M(3) into SO(3)
   It is expecting a 3-by-3 matrix */

  Eigen::Matrix3d Msq = M.transpose() * M;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(Msq, Eigen::ComputeEigenvectors);
  Eigen::Matrix3d u = eigenSolver.eigenvectors().rowwise().reverse();

  int sign = 1;
  if (M.determinant() < 0)
    sign = -1;

  Eigen::Matrix3d dMat = Eigen::Matrix3d::Zero();
  dMat(0, 0) = std::pow(eigenSolver.eigenvalues()(2), -0.5);
  dMat(1, 1) = std::pow(eigenSolver.eigenvalues()(1), -0.5);
  dMat(2, 2) = sign * std::pow(eigenSolver.eigenvalues()(0), -0.5);

  return M * u * dMat * u.transpose();
}

Eigen::Matrix3d meanSO3C(const Eigen::MatrixXd &Rs)
{
  /* Compute the projected mean for a sample of n rotations, Rs.
   This function expects Rs to be a n-by-9 matrix where each row
   represents an observations in SO(3) */

  Eigen::VectorXd Rbarels = Rs.colwise().mean();
  Eigen::Map<Eigen::MatrixXd> Rbar(Rbarels.data(), 3, 3);

  return projectSO3C(Rbar);
}

Eigen::Vector4d exp_quat(const Eigen::Map<Eigen::VectorXd> &x)
{
  Eigen::Vector4d xx(x);
  Eigen::QuaternionMapAlignedd q(xx.data());
  q = expq<double>(q);
  return xx;
}

Eigen::Vector4d log_quat(const Eigen::Map<Eigen::VectorXd> &x)
{
  Eigen::Vector4d xx(x);
  Eigen::QuaternionMapAlignedd q(xx.data());
  q = logq<double>(q);
  return xx;
}

//' @export
// [[Rcpp::export]]
double geodist(const Eigen::Map<Eigen::VectorXd> &x1, const Eigen::Map<Eigen::VectorXd> &x2)
{
  Eigen::Quaterniond q1, q2;
  q1.w() = x1(0);
  q1.x() = x1(1);
  q1.y() = x1(2);
  q1.z() = x1(3);
  q2.w() = x2(0);
  q2.x() = x2(1);
  q2.y() = x2(2);
  q2.z() = x2(3);

  return q1.angularDistance(q2);
}
