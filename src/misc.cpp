#include "misc.h"

double inner_product_with_yinit(const Eigen::Map<Eigen::VectorXd> &q,
                                const Eigen::Map<Eigen::VectorXd> &qinit)
{
  Eigen::Quaterniond qq(q(0), q(1), q(2), q(3));
  Eigen::Quaterniond qqinit(qinit(0), qinit(1), qinit(2), qinit(3));
  return qq._transformVector(Eigen::Vector3d::UnitX()).dot(qqinit._transformVector(Eigen::Vector3d::UnitY()));
}

Rcpp::DataFrame calibrate_xy(const Rcpp::DataFrame &qts,
                             const Eigen::Map<Eigen::VectorXd> &q0)
{
  Eigen::Quaterniond qq0(q0(0), q0(1), q0(2), q0(3));
  Eigen::Vector3d y0 = qq0._transformVector(Eigen::Vector3d::UnitY());
  double inplaneNorm = std::sqrt(y0(0) * y0(0) + y0(1) * y0(1));
  double thetaValue = std::acos(y0(1) / inplaneNorm);

  qq0.w() = std::cos(thetaValue / 2.0);
  qq0.vec() = -std::sin(thetaValue / 2.0) * Eigen::Vector3d::UnitZ();
  y0 = qq0._transformVector(y0);

  if (std::abs(y0(0)) > std::sqrt(std::numeric_limits<double>::epsilon()))
    qq0.vec() *= -1.0;

  unsigned int nGrid = qts.nrows();
  Rcpp::DataFrame outValue = Rcpp::clone(qts);
  Rcpp::NumericVector wValues = outValue["w"];
  Rcpp::NumericVector xValues = outValue["x"];
  Rcpp::NumericVector yValues = outValue["y"];
  Rcpp::NumericVector zValues = outValue["z"];
  Eigen::Quaterniond currentQValue;

  for (unsigned int i = 0;i < nGrid;++i)
  {
    currentQValue.w() = wValues(i);
    currentQValue.x() = xValues(i);
    currentQValue.y() = yValues(i);
    currentQValue.z() = zValues(i);
    currentQValue = qq0 * currentQValue;
    wValues(i) = currentQValue.w();
    xValues(i) = currentQValue.x();
    yValues(i) = currentQValue.y();
    zValues(i) = currentQValue.z();
  }

  return outValue;
}

Eigen::Vector4d rot_q(const Eigen::Map<Eigen::VectorXd> &axis1,
                      const Eigen::Map<Eigen::VectorXd> &axis2)
{
  Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(axis1, axis2);
  return {q.w(), q.x(), q.y(), q.z()};
}
