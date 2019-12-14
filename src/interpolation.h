#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <Rcpp.h>

// [[Rcpp::interfaces(r, cpp)]]

void Normalize(const Rcpp::NumericVector &input,
               Rcpp::NumericVector &output);

Rcpp::NumericVector slerp(const Rcpp::NumericVector &v0,
                          const Rcpp::NumericVector &v1,
                          const double t);

// [[Rcpp::export]]
Rcpp::NumericMatrix RegularizeGrid(const Rcpp::NumericVector &x,
                                   const Rcpp::NumericMatrix &y,
                                   const double xmin,
                                   const double xmax,
                                   const unsigned int outSize = 0);

#endif /* INTERPOLATION_H */
