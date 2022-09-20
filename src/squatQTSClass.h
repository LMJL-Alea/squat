#ifndef SQUATQTSCLASS_H
#define SQUATQTSCLASS_H

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::DataFrame reorient_qts_impl(const Rcpp::DataFrame &qts,
                                  const bool disable_normalization = false);
// [[Rcpp::export]]
Rcpp::DataFrame normalize_qts_impl(const Rcpp::DataFrame &qts);
// [[Rcpp::export]]
Rcpp::DataFrame derivative_qts_impl(const Rcpp::DataFrame &qts);
// [[Rcpp::export]]
Rcpp::DataFrame log_qts_impl(const Rcpp::DataFrame &qts);
// [[Rcpp::export]]
Rcpp::DataFrame exp_qts_impl(const Rcpp::DataFrame &qts);
// [[Rcpp::export]]
Rcpp::List centring_qts_impl(const Rcpp::DataFrame &qts, const bool standardize = false);

#endif /* SQUATQTSCLASS_H */
