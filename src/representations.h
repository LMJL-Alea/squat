#ifndef REPRESENTATIONS_H
#define REPRESENTATIONS_H

#include <RcppEigen.h>

Eigen::MatrixXd GetRotationsFromQuaternions(
    const std::vector<Eigen::VectorXd> &quaternionSample
);

//' @export
// [[Rcpp::export]]
Eigen::Vector4d gmean(
    const std::vector<Eigen::VectorXd> &quaternionSample,
    unsigned int maxIterations = 2000,
    double maxEpsilon = 1.0e-5
);

//' @export
// [[Rcpp::export]]
Eigen::Vector4d gmedian(
    const std::vector<Eigen::VectorXd> &quaternionSample,
    unsigned int maxIterations = 2000,
    double maxEpsilon = 1.0e-5
);

#endif /* REPRESENTATIONS_H */
