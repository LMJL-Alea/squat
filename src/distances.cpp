#include "distances.h"

#include <RcppEigen.h>

Rcpp::NumericMatrix GetCostMatrix(const Rcpp::DataFrame &qts1,
                                  const Rcpp::DataFrame &qts2,
                                  const bool disable_normalization)
{
  unsigned int n1 = qts1.nrows();
  unsigned int n2 = qts2.nrows();
  Rcpp::NumericMatrix costMatrix(n1, n2);
  Eigen::Quaterniond q1Value, q2Value;

  Rcpp::NumericVector w1Values = qts1["w"];
  Rcpp::NumericVector x1Values = qts1["x"];
  Rcpp::NumericVector y1Values = qts1["y"];
  Rcpp::NumericVector z1Values = qts1["z"];

  Rcpp::NumericVector w2Values = qts2["w"];
  Rcpp::NumericVector x2Values = qts2["x"];
  Rcpp::NumericVector y2Values = qts2["y"];
  Rcpp::NumericVector z2Values = qts2["z"];

  for (unsigned int i = 0;i < n1;++i)
  {
    q1Value.w() = w1Values(i);
    q1Value.x() = x1Values(i);
    q1Value.y() = y1Values(i);
    q1Value.z() = z1Values(i);
    if (!disable_normalization)
    {
      q1Value.normalize();
      w1Values(i) = q1Value.w();
      x1Values(i) = q1Value.x();
      y1Values(i) = q1Value.y();
      z1Values(i) = q1Value.z();
    }

    for (unsigned int j = 0;j < n2;++j)
    {
      q2Value.w() = w2Values(j);
      q2Value.x() = x2Values(j);
      q2Value.y() = y2Values(j);
      q2Value.z() = z2Values(j);
      if (!disable_normalization && i == 0)
      {
        q2Value.normalize();
        w2Values(j) = q2Value.w();
        x2Values(j) = q2Value.x();
        y2Values(j) = q2Value.y();
        z2Values(j) = q2Value.z();
      }

      costMatrix(i, j) = q1Value.angularDistance(q2Value);
    }
  }

  return costMatrix;
}
