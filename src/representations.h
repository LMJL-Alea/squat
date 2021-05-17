#ifndef REPRESENTATIONS_H
#define REPRESENTATIONS_H

#include <RcppEigen.h>

// [[Rcpp::interfaces(r, cpp)]]

Eigen::MatrixXd GetRotationsFromQuaternions(
    const std::vector<Eigen::VectorXd> &quaternionSample
);

// [[Rcpp::export]]
Eigen::VectorXd gmean(
    const std::vector<Eigen::VectorXd> &quaternionSample,
    unsigned int maxIterations = 2000,
    double maxEpsilon = 1.0e-5
);

// [[Rcpp::export]]
Eigen::VectorXd gmedian(
    const std::vector<Eigen::VectorXd> &quaternionSample,
    unsigned int maxIterations = 2000,
    double maxEpsilon = 1.0e-5
);

#endif /* REPRESENTATIONS_H */
