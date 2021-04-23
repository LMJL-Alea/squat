#include "representations.h"
#include "rotations.h"

Eigen::MatrixXd GetRotationFromQuaternion(const std::vector<Eigen::VectorXd> &quaternionSample)
{
  unsigned int n = quaternionSample.size();
  Eigen::MatrixXd rotationSample(n, 9);
  Eigen::Quaterniond q;

  for (unsigned int i = 0;i < n;++i)
  {
    q.coeffs() = quaternionSample[i];
    rotationSample.row(i) = Eigen::Map<Eigen::VectorXd>(q.toRotationMatrix().data(), 9);
  }

  return rotationSample;
}

Eigen::VectorXd GetQuaternionFromRotation(const Eigen::Matrix3d &x)
{
  Eigen::Quaterniond q(x);
  return q.coeffs();
}

Eigen::VectorXd GetGeodesicMean(const std::vector<Eigen::VectorXd> &quaternionSample,
                                unsigned int maxIterations,
                                double maxEpsilon)
{
  Eigen::MatrixXd rotationSample = GetRotationFromQuaternion(quaternionSample);
  unsigned int numPoints = rotationSample.rows();
  unsigned int iterations = 0;
  Eigen::MatrixXd Rsi(3, 3);
  Eigen::MatrixXd r;
  Eigen::MatrixXd S = meanSO3C(rotationSample);
  double epsilon = 1.0;

  while (epsilon > maxEpsilon && iterations < maxIterations)
  {
    r = Eigen::MatrixXd::Zero(3, 3);

    for (unsigned int i = 0;i < numPoints;++i)
    {
      for (unsigned int j = 0;j < 9;++j)
        Rsi(j) = rotationSample(i,j);

      r += logSO3C(S.transpose() * Rsi);
    }

    r = r / numPoints;
    S = S * expskewC(r);
    epsilon = r.norm();
    ++iterations;
  }

  return GetQuaternionFromRotation(S);
}
