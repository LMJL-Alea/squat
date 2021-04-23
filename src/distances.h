#ifndef DISTANCES_H
#define DISTANCES_H

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix GetCostMatrix(
    const Rcpp::DataFrame &qts1,
    const Rcpp::DataFrame &qts2,
    const bool disable_normalization = false
);

#endif /* DISTANCES_H */
