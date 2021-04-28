#ifndef MISC_H
#define MISC_H

#include <RcppEigen.h>

//' @export
// [[Rcpp::export]]
double inner_product_with_yinit(
    const Eigen::Map<Eigen::VectorXd> &q,
    const Eigen::Map<Eigen::VectorXd> &qinit
);

//' @export
// [[Rcpp::export]]
Rcpp::DataFrame calibrate_xy(
    const Rcpp::DataFrame &qts,
    const Eigen::Map<Eigen::VectorXd> &q0
);

//' @export
// [[Rcpp::export]]
Eigen::Vector4d rot_q(
    const Eigen::Map<Eigen::VectorXd> &axis1,
    const Eigen::Map<Eigen::VectorXd> &axis2
);

#endif /* MISC_H */
