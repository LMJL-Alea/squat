#ifndef ROTATIONS_H
#define ROTATIONS_H

#include <RcppEigen.h>

Eigen::Matrix3d logSO3C(const Eigen::MatrixXd &R);
Eigen::Matrix3d expskewC(const Eigen::MatrixXd &M);
Eigen::Matrix3d projectSO3C(const Eigen::MatrixXd &M);
Eigen::Matrix3d meanSO3C(const Eigen::MatrixXd &Rs);

#endif /* ROTATIONS_H */
