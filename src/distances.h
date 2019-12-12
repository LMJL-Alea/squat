#pragma once

#include <Rcpp.h>

double GeodesicQuaternionDistance(const Rcpp::NumericMatrix &x,
                                  const Rcpp::NumericMatrix &y,
                                  const unsigned int xIndex,
                                  const unsigned int yIndex);

// [[Rcpp::export]]
Rcpp::NumericMatrix GetCostMatrix(const Rcpp::NumericMatrix &x,
                                  const Rcpp::NumericMatrix &y);

//' L2 Distance for Quaternion Time Series
//'
//' This function computes the L2 distance between two quaternion
//' time series, using the geodesic distance as inter-point
//' quaternion distance.
//'
//' The function currently assumes that the two QTS are evaluated
//' on the same grid.
//'
//' @param x A \code{4 x p} matrix representing the first
//'   quaternion time series.
//' @param y A \code{4 x p} matrix representing the second
//'   quaternion time series.
//'
//' @return A positive scalar providing a measure of distance
//'   between the two input quaternion time series.
//'
//' @export
// [[Rcpp::export]]
double GetL2Distance(const Rcpp::NumericMatrix &x,
                     const Rcpp::NumericMatrix &y);
