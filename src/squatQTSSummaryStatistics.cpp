#include "squatQTSSummaryStatistics.h"
#include "rotations.h"

Eigen::MatrixXd GetRotationsFromQuaternions(const std::vector<Eigen::VectorXd> &quaternionSample)
{
  unsigned int nQuaternions = quaternionSample.size();
  Eigen::MatrixXd rotationSample(nQuaternions, 9);
  Eigen::Quaterniond qValue;
  Eigen::Vector4d workValue;

  for (unsigned int i = 0;i < nQuaternions;++i)
  {
    workValue = quaternionSample[i];
    qValue = Eigen::Quaterniond(workValue(0), workValue(1), workValue(2), workValue(3));
    rotationSample.row(i) = Eigen::Map<Eigen::VectorXd>(qValue.toRotationMatrix().data(), 9);
  }

  return rotationSample;
}

Eigen::VectorXd gmean(const std::vector<Eigen::VectorXd> &quaternionSample,
                      unsigned int maxIterations,
                      double maxEpsilon)
{
  Eigen::MatrixXd rotationSample = GetRotationsFromQuaternions(quaternionSample);
  unsigned int numPoints = rotationSample.rows();
  unsigned int iterations = 0;
  Eigen::Matrix3d Rsi;
  Eigen::Matrix3d r;
  Eigen::Matrix3d S = meanSO3C(rotationSample);
  double epsilon = 1.0;

  while (epsilon > maxEpsilon && iterations < maxIterations)
  {
    r = Eigen::Matrix3d::Zero();

    for (unsigned int i = 0;i < numPoints;++i)
    {
      for (unsigned int j = 0;j < 9;++j)
        Rsi(j) = rotationSample(i, j);

      r += logSO3C(S.transpose() * Rsi);
    }

    r = r / numPoints;
    S = S * expskewC(r);
    epsilon = r.norm();
    ++iterations;
  }

  Eigen::Quaterniond meanQValue(S);
  Eigen::Vector4d outValue;
  outValue(0) = meanQValue.w();
  outValue(1) = meanQValue.x();
  outValue(2) = meanQValue.y();
  outValue(3) = meanQValue.z();
  return outValue;
}

double gvariance(const std::vector<Eigen::VectorXd> &quaternionSample,
                 const Eigen::VectorXd &quaternionMean)
{
  unsigned int numPoints = quaternionSample.size();
  Eigen::Quaterniond workQuaternion, meanQuaternion(quaternionMean(0), quaternionMean(1), quaternionMean(2), quaternionMean(3));

  double varValue = 0.0;
  for (unsigned int i = 0;i < numPoints;++i)
  {
    workQuaternion = Eigen::Quaterniond(quaternionSample[i](0), quaternionSample[i](1), quaternionSample[i](2), quaternionSample[i](3));
    double distValue = workQuaternion.angularDistance(meanQuaternion);
    varValue += distValue * distValue;
  }

  return varValue;
}

Eigen::VectorXd gmedian(const std::vector<Eigen::VectorXd> &quaternionSample,
                        const unsigned int maxIterations,
                        const double maxEpsilon)
{
  Eigen::MatrixXd rotationSample = GetRotationsFromQuaternions(quaternionSample);
  unsigned int numPoints = rotationSample.rows();
  unsigned int iterations = 0;
  Eigen::Matrix3d S = meanSO3C(rotationSample);
  Eigen::Matrix3d Snew;
  Eigen::Matrix3d delta;
  Eigen::Matrix3d Rsi;
  Eigen::Matrix3d vi;
  double epsilon = 1.0;

  while (epsilon > maxEpsilon && iterations < maxIterations)
  {
    delta = Eigen::Matrix3d::Zero();
    double denom = 0;

    for (unsigned int i = 0;i < numPoints;++i)
    {
      for (unsigned int j = 0;j < 9;++j)
        Rsi(j) = rotationSample(i, j);

      vi = logSO3C(Rsi * S.transpose());
      double vin = std::max(vi.operatorNorm(), 1.0e-5);

      vin = 1.0 / vin;
      delta += vi * vin;
      denom += vin;
    }

    delta /= denom;
    Snew = expskewC(delta) * S;

    ++iterations;
    epsilon = (Snew - S).operatorNorm();
    S = Snew;
  }

  Eigen::Quaterniond medianQValue(S);
  Eigen::Vector4d outValue;
  outValue(0) = medianQValue.w();
  outValue(1) = medianQValue.x();
  outValue(2) = medianQValue.y();
  outValue(3) = medianQValue.z();
  return outValue;
}
