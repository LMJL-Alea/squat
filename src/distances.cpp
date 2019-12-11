#include <Rcpp.h>

//' @export
// [[Rcpp::export]]
double GetSquaredL2Distance(const Rcpp::NumericMatrix &x, const Rcpp::NumericMatrix &y)
{
  unsigned int nx = x.ncol();
  unsigned int ny = y.ncol();

  if (nx != ny)
    Rcpp::stop("Evaluation of the L2 Distance requires curves to be evaluated on the same grid.");

  double sqDistance = 0.0;
  for (unsigned int i = 0;i < nx;++i)
  {
    double realPart = 0.0;
    for (unsigned int j = 0;j < 4;++j)
      realPart += x(j, i) * y(j, i);
    realPart = std::abs(realPart);

    const double DOT_THRESHOLD = 0.9995;
    if (realPart > DOT_THRESHOLD)
      continue;

    sqDistance += 2.0 * std::acos(realPart);
  }

  return sqDistance;
}
