#include "distances.h"

double GeodesicQuaternionDistance(const Rcpp::NumericMatrix &x,
                                  const Rcpp::NumericMatrix &y,
                                  const unsigned int xIndex,
                                  const unsigned int yIndex)
{
  double realPart = 0.0;

  for (unsigned int j = 0;j < 4;++j)
    realPart += x(j, xIndex) * y(j, yIndex);

  realPart = std::abs(realPart);

  const double DOT_THRESHOLD = 0.9995;

  if (realPart > DOT_THRESHOLD)
    return 0.0;

  return 2.0 * std::acos(realPart);
}

Rcpp::NumericMatrix GetCostMatrix(const Rcpp::NumericMatrix &x,
                                  const Rcpp::NumericMatrix &y)
{
  unsigned int nx = x.ncol();
  unsigned int ny = y.ncol();
  Rcpp::NumericMatrix costMatrix(nx, ny);

  for (unsigned int i = 0;i < nx;++i)
    for (unsigned int j = 0;j < ny;++j)
      costMatrix(i, j) = GeodesicQuaternionDistance(x, y, i, j);

  return costMatrix;
}

double GetL2Distance(const Rcpp::NumericMatrix &x,
                     const Rcpp::NumericMatrix &y)
{
  unsigned int nx = x.ncol();
  unsigned int ny = y.ncol();

  if (nx != ny)
    Rcpp::stop("Evaluation of the L2 Distance requires curves to be evaluated on the same grid.");

  double sqDistance = 0.0;
  for (unsigned int i = 0;i < nx;++i)
  {
    double distValue = GeodesicQuaternionDistance(x, y, i, i);
    sqDistance += distValue * distValue;
  }

  return std::sqrt(sqDistance);
}
