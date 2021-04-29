#include "rotations.h"

Eigen::Matrix3d logSO3C(const Eigen::Matrix3d &R)
{
  /* This function takes a 3-by-3 rotation matrix (in SO(3)) and
   returns its logarithm, a 3-by-3 skew symmetric matrix (in so(3)) */

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

Eigen::Matrix3d expskewC(const Eigen::Matrix3d &M)
{
  /* This function takes a 3-by-3 skew symmetric matrix (in so(3)) and
   returns its exponential, a 3-by-3 rotation matrix (in SO(3)) */

  Eigen::Matrix3d expM = Eigen::Matrix3d::Identity();

  double a = std::sqrt(0.5 * (M.transpose() * M).trace());

  if (std::abs(a) < std::sqrt(std::numeric_limits<double>::epsilon()))
    return expM;

  expM += (std::sin(a) / a) * M + (1.0 - std::cos(a)) * M * M / (a * a);

  return expM;
}

Eigen::Matrix3d projectSO3C(const Eigen::Matrix3d &M)
{
  /* This function projects a 3-by-3 symmetric matrix (in M(3)) into SO(3). */

  Eigen::Matrix3d Msq = M.transpose() * M;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(Msq, Eigen::ComputeEigenvectors);
  Eigen::Matrix3d u = eigenSolver.eigenvectors().rowwise().reverse();

  Eigen::Matrix3d dMat = eigenSolver.eigenvalues().array().rsqrt().matrix().reverse().asDiagonal();
  if (M.determinant() < 0)
    dMat(2, 2) *= -1.0;

  return M * u * dMat * u.transpose();
}

Eigen::Matrix3d meanSO3C(const Eigen::MatrixXd &Rs)
{
  /* Compute the projected mean for a sample of n rotations, Rs.
   This function expects Rs to be a n-by-9 matrix where each row
   represents an observation in SO(3) */

  Eigen::VectorXd Rbarels = Rs.colwise().mean();
  Eigen::Map<Eigen::Matrix3d> Rbar(Rbarels.data());

  return projectSO3C(Rbar);
}
