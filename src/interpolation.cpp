#include <Rcpp.h>

void Normalize(const Rcpp::NumericVector &input, Rcpp::NumericVector &output)
{
  double norm = 0.0;
  for (unsigned int i = 0;i < 4;++i)
    norm += input[i] * input[i];
  norm = std::sqrt(norm);
  output = input / norm;
}

// [[Rcpp::export]]
Rcpp::NumericVector slerp(const Rcpp::NumericVector &v0, const Rcpp::NumericVector &v1, const double t)
{
  // Only unit quaternions are valid rotations.
  // Normalize to avoid undefined behavior.
  Rcpp::NumericVector q0, q1;
  Normalize(v0, q0);
  Normalize(v1, q1);

  // Compute the cosine of the angle between the two vectors.
  double dot = 0.0;
  for (unsigned int i = 0;i < 4;++i)
    dot += q0[i] * q1[i];

  // If the dot product is negative, slerp won't take
  // the shorter path. Note that v1 and -v1 are equivalent when
  // the negation is applied to all four components. Fix by
  // reversing one quaternion.
  if (dot < 0.0)
  {
    q1 = q1 * (-1.0);
    dot = -dot;
  }

  const double DOT_THRESHOLD = 0.9995;
  if (dot > DOT_THRESHOLD)
  {
    // If the inputs are too close for comfort, linearly interpolate
    // and normalize the result.
    Rcpp::NumericVector result = q0 + t * (q1 - q0);
    Normalize(result, result);
    return result;
  }

  // Since dot is in range [0, DOT_THRESHOLD], acos is safe
  double theta_0 = std::acos(dot);        // theta_0 = angle between input vectors
  double theta = theta_0 * t;             // theta = angle between v0 and result
  double sin_theta = std::sin(theta);     // compute this value only once
  double sin_theta_0 = std::sin(theta_0); // compute this value only once

  double s0 = std::cos(theta) - dot * sin_theta / sin_theta_0;  // == sin(theta_0 - theta) / sin(theta_0)
  double s1 = sin_theta / sin_theta_0;

  return (s0 * q0) + (s1 * q1);
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix RegularizeGrid(const Rcpp::NumericVector &x, const Rcpp::NumericMatrix &y, const unsigned int outSize = 0)
{
  // Assumes x is sorted in ascending order
  unsigned int sizeIn = x.size();
  unsigned int sizeOut = (outSize == 0) ? sizeIn : outSize;
  unsigned int posInf = 0;
  unsigned int posSup = sizeIn - 1;
  double xmin = x[posInf];
  double xmax = x[posSup];
  double step = (xmax - xmin) / (sizeOut - 1.0);
  Rcpp::NumericVector Qinf, Qsup;
  Rcpp::NumericMatrix yOut(4, sizeOut);

  for (unsigned int i = 0;i < sizeOut;++i)
  {
    double newx = xmin + (double)i * step;

    double xinf = x[posInf];
    while (xinf <= newx & posInf < sizeIn)
    {
      posInf += 1;
      xinf = x[posInf];
    }
    if (posInf > 0)
      posInf -= 1;
    if (posInf >= sizeIn)
      posInf = sizeIn - 1;
    Qinf = y(Rcpp::_, posInf);

    double xsup = x[posSup];
    while (xsup > newx & posSup >= 0)
    {
      posSup -= 1;
      xsup = x[posSup];
    }
    if (posSup < sizeIn - 1)
      posSup += 1;
    if (posSup < 0)
      posSup = 0;
    Qsup = y(Rcpp::_, posSup);

    if (xsup == xinf)
      yOut(Rcpp::_, i) = Qinf;
    else
    {
      double p = (xsup - newx) / (xsup - xinf);
      yOut(Rcpp::_, i) = slerp(Qinf, Qsup, p);
    }
  }

  return yOut;
}
