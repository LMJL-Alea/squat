#ifndef ROTATIONS_HPP
#define ROTATIONS_HPP

#include "rotations.h"

template<typename T>
Eigen::Quaternion<T> expq(const Eigen::Quaternion<T>& q)
{
  T a = q.vec().norm();
  T exp_w = std::exp(q.w());

  if (a == T(0))
    return Eigen::Quaternion<T>(exp_w, T(0), T(0), T(0));

  Eigen::Quaternion<T> res;
  res.w() = exp_w * T(std::cos(a));
  res.vec() = exp_w * T(std::sin(a) / a) * q.vec();

  return res;
}

template<typename T>
Eigen::Quaternion<T> logq(const Eigen::Quaternion<T>& q)
{
  T exp_w = q.norm();
  T w = std::log(exp_w);
  T a = std::acos(q.w() / exp_w);

  if (a == T(0))
  {
    return Eigen::Quaternion<T>(w, T(0), T(0), T(0));
  }

  Eigen::Quaternion<T> res;
  res.w() = w;
  res.vec() = q.vec() / exp_w / (sin(a) / a);

  return res;
}

#endif /* ROTATIONS_HPP */
