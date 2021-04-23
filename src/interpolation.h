#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <Rcpp.h>

// [[Rcpp::interfaces(r, cpp)]]

//' @export
// [[Rcpp::export]]
Rcpp::DataFrame resample_qts(
    const Rcpp::DataFrame &qts,
    const unsigned int nout = 0,
    const bool disable_normalization = false
);

#endif /* INTERPOLATION_H */
