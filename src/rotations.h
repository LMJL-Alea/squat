#ifndef ROTATIONS_H
#define ROTATIONS_H

#include <RcppEigen.h>

Eigen::Matrix3d logSO3C(const Eigen::Matrix3d &R);
Eigen::Matrix3d expskewC(const Eigen::Matrix3d &M);
Eigen::Matrix3d projectSO3C(const Eigen::Matrix3d &M);
Eigen::Matrix3d meanSO3C(const Eigen::MatrixXd &Rs);

Eigen::Quaterniond expq(const Eigen::Quaterniond& q);
Eigen::Quaterniond logq(const Eigen::Quaterniond& q);

// [[Rcpp::export]]
Eigen::Vector4d geometric_mean(
    const std::vector<Eigen::VectorXd> &quaternionSample,
    unsigned int maxIterations = 2000,
    double maxEpsilon = 1.0e-5
);

#endif /* ROTATIONS_H */
