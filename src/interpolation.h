#pragma once

#include <Rcpp.h>

void Normalize(const Rcpp::NumericVector &input,
               Rcpp::NumericVector &output);

Rcpp::NumericVector slerp(const Rcpp::NumericVector &v0,
                          const Rcpp::NumericVector &v1,
                          const double t);

//' Grid regularization
//'
//' This function makes sure that a quaternion time series
//' is evaluated on a grid with fixed step size.
//'
//' @param x A numeric vector providing the original evaluation grid.
//' @param y A \code{4 x length(x)} matrix providing the original QTS.
//' @param outSize An integer specifying the size of the output
//'   evaluation grid. By default, it takes the same size as the input
//'   evaluation grid.
//'
//' @return A \code{4 x length(x)} matrix providing the regularized QTS.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix RegularizeGrid(const Rcpp::NumericVector &x,
                                   const Rcpp::NumericMatrix &y,
                                   const unsigned int outSize = 0);
