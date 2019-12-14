#ifndef REPRESENTATIONS_H
#define REPRESENTATIONS_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::interfaces(r, cpp)]]

Rcpp::NumericMatrix GetAxisFromQuaternion(const Rcpp::NumericMatrix &quaternionSample);
Rcpp::NumericMatrix GetAxisFromRotation(const Rcpp::NumericMatrix &rotationSample);

Rcpp::NumericVector GetAngleFromQuaternion(const Rcpp::NumericMatrix &quaternionSample);
Rcpp::NumericVector GetAngleFromRotation(const Rcpp::NumericMatrix &rotationSample);

Rcpp::NumericMatrix GetRotationFromQuaternion(const Rcpp::NumericMatrix &quaternionSample);
Rcpp::NumericMatrix GetQuaternionFromRotation(const Rcpp::NumericMatrix &rotationSample);

// [[Rcpp::export]]
Rcpp::NumericMatrix GetGeodesicMean(const Rcpp::NumericMatrix &quaternionSample,
                                    unsigned int maxIterations = 2000,
                                    double maxEpsilon = 1.0e-5);

#endif /* REPRESENTATIONS_H */
