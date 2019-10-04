#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix GetCostMatrix(const Rcpp::NumericMatrix &x, const Rcpp::NumericMatrix &y) {
  unsigned int nx = x.ncol();
  unsigned int ny = y.ncol();
  Rcpp::NumericMatrix costMatrix(nx, ny);
  for (unsigned int i = 0;i < nx;++i)
  {
    for (unsigned int j = 0;j < ny;++j)
    {
      double realPart = 0.0;
      for (unsigned int k = 0;k < 4;++k)
        realPart += x(k, i) * y(k, j);
      realPart = std::abs(realPart);
      if (realPart >  1.0)
        realPart =  1.0;
      costMatrix(i, j) = 2.0 * std::acos(realPart);
    }
  }
  return costMatrix;
}
