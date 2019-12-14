// [[Rcpp::depends(rotations)]]
#include <rotations.h>

#include "representations.h"

Rcpp::NumericMatrix GetAxisFromQuaternion(const Rcpp::NumericMatrix &quaternionSample)
{
  Rcpp::NumericVector angleValues = GetAngleFromQuaternion(quaternionSample);

  unsigned int numPoints = quaternionSample.nrow();
  Rcpp::NumericMatrix axes(numPoints, 3);

  for (unsigned int i = 0;i < numPoints;++i)
  {
    double angle = angleValues[i];
    double sinValue = std::sin(angle / 2.0);

    for (unsigned int j = 0;j < 3;++j)
    {
      if (std::abs(sinValue) < std::numeric_limits<double>::epsilon())
        axes(i, j) = 0.0;
      else
        axes(i, j) = quaternionSample(i, j + 1) / sinValue;
    }
  }

  return axes;
}

Rcpp::NumericMatrix GetAxisFromRotation(const Rcpp::NumericMatrix &rotationSample)
{
  unsigned int numPoints = rotationSample.nrow();
  Rcpp::NumericMatrix axes(numPoints, 3);
  Rcpp::NumericMatrix workRotation(3, 3);

  for (unsigned int i = 0;i < numPoints;++i)
  {
    double trace = rotationSample(i, 0) + rotationSample(i, 4) + rotationSample(i, 8);

    if (std::abs(3.0 - trace) < std::numeric_limits<double>::epsilon())
      continue;

    int pos = 0;

    for (unsigned int j = 0;j < 3;++j)
    {
      for (unsigned int k = 0;k < 3;++k)
      {
        workRotation(k, j) = rotationSample(i, pos);
        ++pos;
      }
    }

    pos = 2;
    double signValue = -1.0;
    double sqNorm = 0.0;
    for (unsigned int j = 0;j < 3;++j)
    {
      for (unsigned int k = 0;k < 3;++k)
      {
        if (j >= k)
          continue;

        double workValue = signValue * (workRotation(j, k) - workRotation(k, j));
        axes(i, pos) = workValue;
        signValue *= -1.0;
        sqNorm += workValue * workValue;
        --pos;
      }
    }

    if (sqNorm > std::numeric_limits<double>::epsilon())
    {
      for (unsigned int j = 0;j < 3;++j)
        axes(i, j) /= std::sqrt(sqNorm);
    }
  }

  return axes;
}

Rcpp::NumericVector GetAngleFromQuaternion(const Rcpp::NumericMatrix &quaternionSample)
{
  unsigned int numPoints = quaternionSample.nrow();
  Rcpp::NumericVector angles(numPoints);

  for (unsigned int i = 0;i < numPoints;++i)
    angles[i] = 2.0 * std::acos(quaternionSample(i, 0));

  return angles;
}

Rcpp::NumericVector GetAngleFromRotation(const Rcpp::NumericMatrix &rotationSample)
{
  return Rcpp::wrap(rotations::rdistSO3C(
    Rcpp::as<arma::mat>(rotationSample),
    Rcpp::as<arma::mat>(Rcpp::NumericMatrix::diag(3, 1)))
  );
}

Rcpp::NumericMatrix GetRotationFromQuaternion(const Rcpp::NumericMatrix &quaternionSample)
{
  Rcpp::NumericVector angles = GetAngleFromQuaternion(quaternionSample);
  Rcpp::NumericMatrix axes = GetAxisFromQuaternion(quaternionSample);

  unsigned int numPoints = quaternionSample.nrow();
  Rcpp::NumericMatrix rotationSample(numPoints, 9);
  arma::mat workMatrix;
  arma::rowvec workVector;

  for (unsigned int i = 0;i < numPoints;++i)
  {
    workVector = axes.row(i);
    if (arma::norm(workVector) < std::numeric_limits<double>::epsilon())
    {
      rotationSample(i, 0) = 1.0;
      rotationSample(i, 4) = 1.0;
      rotationSample(i, 8) = 1.0;
      continue;
    }

    workMatrix = rotations::eskewC(workVector);
    double thetaValue = angles[i];
    double cosValue = std::cos(thetaValue);
    double sinValue = std::sin(thetaValue);
    unsigned int pos = 0;

    for (unsigned int k = 0;k < 3;++k)
    {
      for (unsigned int j = 0;j < 3;++j)
      {
        double projectionValue = axes(i, j) * axes(i, k);
        double workValue = projectionValue;
        workValue += ((double)(j == k) - projectionValue) * cosValue;
        workValue += workMatrix(j, k) * sinValue;
        rotationSample(i, pos) = workValue;
        ++pos;
      }
    }
  }

  return rotationSample;
}

Rcpp::NumericMatrix GetQuaternionFromRotation(const Rcpp::NumericMatrix &rotationSample)
{
  Rcpp::NumericVector angles = GetAngleFromRotation(rotationSample);
  Rcpp::NumericMatrix axes = GetAxisFromRotation(rotationSample);

  unsigned int numPoints = rotationSample.nrow();
  Rcpp::NumericMatrix quaternionSample(numPoints, 4);

  for (unsigned int i = 0;i < numPoints;++i)
  {
    double angleValue = angles[i] / 2.0;
    quaternionSample(i, 0) = std::cos(angleValue);
    double sinValue = std::sin(angleValue);

    for (unsigned int j = 0;j < 3;++j)
      quaternionSample(i, j + 1) = sinValue * axes(i, j);
  }

  return quaternionSample;
}

Rcpp::NumericMatrix GetGeodesicMean(const Rcpp::NumericMatrix &quaternionSample,
                                    unsigned int maxIterations,
                                    double maxEpsilon)
{
  Rcpp::NumericMatrix rotationSample = GetRotationFromQuaternion(quaternionSample);
  unsigned int numPoints = rotationSample.nrow();
  unsigned int iterations = 0;
  arma::mat33 Rsi, r;
  arma::mat S = rotations::meanSO3C(Rcpp::as<arma::mat>(rotationSample));
  double epsilon = 1.0;
  Rsi.zeros();

  while (epsilon > maxEpsilon && iterations < maxIterations)
  {
    r.zeros();

    for (unsigned int i = 0;i < numPoints;++i)
    {
      for (unsigned int j = 0;j < 9;++j)
        Rsi(j) = rotationSample(i,j);

      r = r + rotations::logSO3C(S.t() * Rsi);
    }

    r = r / numPoints;
    S = S * rotations::expskewC(r);
    epsilon = arma::norm(r, "fro");
    ++iterations;
  }

  Rcpp::NumericMatrix rotationMean(1, 9);

  unsigned int pos = 0;
  for (unsigned int j = 0;j < 3;++j)
  {
    for (unsigned int k = 0;k < 3;++k)
    {
      rotationMean(0, pos) = S(k, j);
      ++pos;
    }
  }

  return GetQuaternionFromRotation(rotationMean);
}
