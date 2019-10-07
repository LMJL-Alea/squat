#include <Rcpp.h>

void Normalize(const Rcpp::NumericVector &input, Rcpp::NumericVector &output)
{
  double norm = 0.0;
  for (unsigned int i = 0;i < 4;++i)
    norm += input[i] * input[i];
  norm = std::sqrt(norm);
  output = input / norm;
}

Rcpp::NumericVector slerp(const Rcpp::NumericVector &v0, const Rcpp::NumericVector &v1, const double t) {
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

// [[Rcpp::export]]
Rcpp::NumericMatrix RegularizeGrid(const Rcpp::NumericVector &x, const Rcpp::NumericMatrix &y, const double step = 1)
{
  unsigned int gridSize = x.size();
  double xmin = x[0], xmax = x[gridSize - 1];
  std::vector<Rcpp::NumericVector> interpolatedValues;
  Rcpp::NumericVector Qinf, Qsup;
  double newx = xmin;

  while (newx < xmax)
  {
    unsigned int pos = 0;
    double xinf = x[pos];
    while (xinf < newx)
    {
      pos += 1;
      xinf = x[pos];
    }
    if (pos > 0)
    {
      pos -= 1;
      xinf = x[pos];
    }
    Qinf = y(Rcpp::_, pos);

    pos = gridSize - 1;
    double xsup = x[pos];
    while (xsup > newx)
    {
      pos -= 1;
      xsup = x[pos];
    }
    if (pos < gridSize - 1)
    {
      pos += 1;
      xsup = x[pos];
    }
    Qsup = y(Rcpp::_, pos);

    double p = (xsup - newx) / (xsup - xinf);
    interpolatedValues.push_back(slerp(Qinf, Qsup, p));

    newx += step;
  }

  gridSize = interpolatedValues.size();
  Rcpp::NumericMatrix yreg(4, gridSize);
  for (unsigned int i = 0;i < 4;++i)
    for (unsigned int j = 0;j < gridSize;++j)
      yreg(i, j) = interpolatedValues[j][i];

  return yreg;
}
