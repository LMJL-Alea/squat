#ifndef ROTATIONS_H
#define ROTATIONS_H

#include <RcppEigen.h>

Eigen::Matrix3d logSO3C(const Eigen::MatrixXd &R);
Eigen::Matrix3d expskewC(const Eigen::MatrixXd &M);
Eigen::Matrix3d projectSO3C(const Eigen::MatrixXd &M);
Eigen::Matrix3d meanSO3C(const Eigen::MatrixXd &Rs);

template<typename T>
Eigen::Quaternion<T> expq(const Eigen::Quaternion<T>& q);
template<typename T>
Eigen::Quaternion<T> logq(const Eigen::Quaternion<T>& q);

//' @export
// [[Rcpp::export]]
Eigen::Vector4d exp_quat(const Eigen::Map<Eigen::VectorXd> &x);

//' @export
// [[Rcpp::export]]
Eigen::Vector4d log_quat(const Eigen::Map<Eigen::VectorXd> &x);

#endif /* ROTATIONS_H */

#include "rotations.hpp"
