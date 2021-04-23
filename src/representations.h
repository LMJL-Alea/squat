#ifndef REPRESENTATIONS_H
#define REPRESENTATIONS_H

#include <RcppEigen.h>

// [[Rcpp::interfaces(r, cpp)]]

Eigen::MatrixXd GetRotationFromQuaternion(
    const std::vector<Eigen::VectorXd> &quaternionSample
);

Eigen::VectorXd GetQuaternionFromRotation(const Eigen::Matrix3d &x);

//' @export
// [[Rcpp::export]]
Eigen::VectorXd GetGeodesicMean(
    const std::vector<Eigen::VectorXd> &quaternionSample,
    unsigned int maxIterations = 2000,
    double maxEpsilon = 1.0e-5
);

#endif /* REPRESENTATIONS_H */
