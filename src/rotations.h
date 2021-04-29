#ifndef ROTATIONS_H
#define ROTATIONS_H

#include <RcppEigen.h>

Eigen::Matrix3d logSO3C(const Eigen::Matrix3d &R);
Eigen::Matrix3d expskewC(const Eigen::Matrix3d &M);
Eigen::Matrix3d projectSO3C(const Eigen::Matrix3d &M);
Eigen::Matrix3d meanSO3C(const Eigen::MatrixXd &Rs);

template<typename T>
Eigen::Quaternion<T> expq(const Eigen::Quaternion<T>& q);
template<typename T>
Eigen::Quaternion<T> logq(const Eigen::Quaternion<T>& q);

#endif /* ROTATIONS_H */

#include "rotations.hpp"
