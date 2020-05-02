#include "interpolation.h"

void Normalize(const Rcpp::NumericVector &input,
               Rcpp::NumericVector &output)
{
  double norm = 0.0;
  for (unsigned int i = 0;i < 4;++i)
    norm += input[i] * input[i];
  norm = std::sqrt(norm);
  output = input / norm;
}

Rcpp::NumericVector slerp(const Rcpp::NumericVector &v0,
                          const Rcpp::NumericVector &v1,
                          const double t)
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
    q0 += t * (q1 - q0);
    Normalize(q0, q0);
    return q0;
  }

  // Since dot is in range [0, DOT_THRESHOLD], acos is safe
  double theta_0 = std::acos(dot);        // theta_0 = angle between input vectors
  double theta = theta_0 * t;             // theta = angle between v0 and result
  double sin_theta = std::sin(theta);     // compute this value only once
  double sin_theta_0 = std::sin(theta_0); // compute this value only once

  double s0 = std::cos(theta) - dot * sin_theta / sin_theta_0;  // == sin(theta_0 - theta) / sin(theta_0)
  double s1 = sin_theta / sin_theta_0;

  q0 = (s0 * q0) + (s1 * q1);
  Normalize(q0, q0);
  return q0;
}

// Rcpp::NumericMatrix RegularizeGrid(const Rcpp::NumericVector &x,
//                                    const Rcpp::NumericMatrix &y,
//                                    const double xmin,
//                                    const double xmax,
//                                    const unsigned int outSize)
// {
//   // Assumes x is sorted in ascending order
//   unsigned int sizeIn = x.size();
//   unsigned int sizeOut = (outSize == 0) ? sizeIn : outSize;
//   int posInf = 0;
//   int posSup = sizeIn - 1;
//   double step = (xmax - xmin) / (sizeOut - 1.0);
//   Rcpp::NumericVector Qinf, Qsup;
//   Rcpp::NumericMatrix yOut(4, sizeOut);
//
//   for (unsigned int i = 0;i < sizeOut;++i)
//   {
//     double newx = xmin + (double)i * step;
//
//     posInf = std::lower_bound(x.begin(), x.end(), newx) - x.begin();
//     if (posInf == 0) // first element is already greater or equal
//     {
//       if (std::abs(x[0] - newx) < std::numeric_limits<double>::epsilon())
//         posInf = 0;
//       else
//         Rcpp::stop("inf not allowed");
//     }
//     else if (posInf == sizeIn) // last element is still smaller
//       posInf = sizeIn - 1;
//     else if (x[posInf] > newx)
//       posInf -= 1;
//     double xinf = x[posInf];
//     Qinf = y(Rcpp::_, posInf);
//
//     posSup = std::upper_bound(x.begin(), x.end(), newx) - x.begin();
//     if (posSup == sizeIn) // all elements are smaller or equal
//     {
//       if (std::abs(x[sizeIn - 1] - newx) < std::numeric_limits<double>::epsilon())
//         posSup = sizeIn - 1;
//       else
//         Rcpp::stop("sup not allowed");
//     }
//     else if (posSup == 0) // first element is still larger
//       posSup = 0;
//     else if (x[posSup - 1] == newx)
//       posSup -= 1;
//     double xsup = x[posSup];
//     Qsup = y(Rcpp::_, posSup);
//
//     double range = xsup - xinf;
//
//     if (xsup <= xinf)
//     {
//       Normalize(Qsup, Qsup);
//       yOut(Rcpp::_, i) = Qsup;
//     }
//     else
//     {
//       double p = (newx - xinf) / range;
//       yOut(Rcpp::_, i) = slerp(Qinf, Qsup, p);
//     }
//   }
//
//   return yOut;
// }

Rcpp::NumericMatrix RegularizeGrid(const Rcpp::NumericVector &x,
                                   const Rcpp::NumericMatrix &y,
                                   const double xmin,
                                   const double xmax,
                                   const unsigned int outSize)
{
  // Assumes x is sorted in ascending order
  int sizeIn = x.size();
  unsigned int sizeOut = (outSize == 0) ? sizeIn : outSize;
  int posInf = 0;
  int posSup = sizeIn - 1;
  double step = (xmax - xmin) / (sizeOut - 1.0);
  Rcpp::NumericVector Qinf, Qsup;
  Rcpp::NumericMatrix yOut(4, sizeOut);
  double epsValue = std::sqrt(std::numeric_limits<double>::epsilon());

  for (unsigned int i = 0;i < sizeOut;++i)
  {
    double newx = xmin + (double)i * step;

    if (std::abs(newx - x[0]) < epsValue)
      newx = x[0];

    if (std::abs(newx - x[sizeIn - 1]) < epsValue)
      newx = x[sizeIn - 1];

    if (newx < x[0] || newx > x[sizeIn - 1])
    {
      yOut(Rcpp::_, i) = Rcpp::rep(R_NaN, 4);
      continue;
    }

    double xinf = x[posInf];
    while (xinf <= newx && posInf < sizeIn)
    {
      posInf += 1;
      xinf = x[posInf];
    }
    if (xinf > newx && posInf > 0)
      posInf -= 1;
    else if (posInf >= sizeIn)
      posInf = sizeIn - 1;
    xinf = x[posInf];
    Qinf = y(Rcpp::_, posInf);

    posSup = sizeIn - 1;
    double xsup = x[posSup];
    while (xsup >= newx && posSup >= 0)
    {
      posSup -= 1;
      xsup = x[posSup];
    }
    if (xsup < newx && posSup < sizeIn - 1)
      posSup += 1;
    else if (posSup < 0)
      posSup = 0;
    xsup = x[posSup];
    Qsup = y(Rcpp::_, posSup);

    double range = xsup - xinf;

    if (xsup <= xinf)
    {
      Normalize(Qinf, Qinf);
      yOut(Rcpp::_, i) = Qinf;
    }
    else
    {
      double p = (newx - xinf) / range;
      yOut(Rcpp::_, i) = slerp(Qinf, Qsup, p);
    }
  }

  return yOut;
}
