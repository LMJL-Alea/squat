#include "squatQTSClass.h"
#include "squatSO3Utils.h"
#include <RcppEigen.h>

Rcpp::DataFrame reorient_qts_impl(const Rcpp::DataFrame &qts,
                                  const bool disable_normalization)
{
  unsigned int nSamples = qts.nrows();
  Eigen::Quaterniond qValue;
  Rcpp::DataFrame resValue = Rcpp::clone(qts);
  Rcpp::NumericVector wValues = resValue["w"];
  Rcpp::NumericVector xValues = resValue["x"];
  Rcpp::NumericVector yValues = resValue["y"];
  Rcpp::NumericVector zValues = resValue["z"];

  Eigen::Quaterniond refValue;
  refValue.w() = wValues(0);
  refValue.x() = xValues(0);
  refValue.y() = yValues(0);
  refValue.z() = zValues(0);
  if (!disable_normalization)
    refValue.normalize();
  refValue = refValue.inverse();

  for (unsigned int i = 0;i < nSamples;++i)
  {
    qValue.w() = wValues(i);
    qValue.x() = xValues(i);
    qValue.y() = yValues(i);
    qValue.z() = zValues(i);
    if (!disable_normalization)
      qValue.normalize();

    qValue = refValue * qValue;

    wValues(i) = qValue.w();
    xValues(i) = qValue.x();
    yValues(i) = qValue.y();
    zValues(i) = qValue.z();
  }

  resValue.attr("class") = Rcpp::CharacterVector::create("qts", "tbl_df", "tbl", "data.frame");
  return resValue;
}

Rcpp::DataFrame normalize_qts_impl(const Rcpp::DataFrame &qts)
{
  unsigned int nGrid = qts.nrows();
  Rcpp::DataFrame outValue = Rcpp::clone(qts);
  Rcpp::NumericVector wValues = outValue["w"];
  Rcpp::NumericVector xValues = outValue["x"];
  Rcpp::NumericVector yValues = outValue["y"];
  Rcpp::NumericVector zValues = outValue["z"];
  Eigen::Quaterniond qValue;

  for (unsigned int i = 0;i < nGrid;++i)
  {
    qValue.w() = wValues(i);
    qValue.x() = xValues(i);
    qValue.y() = yValues(i);
    qValue.z() = zValues(i);
    qValue.normalize();
    wValues(i) = qValue.w();
    xValues(i) = qValue.x();
    yValues(i) = qValue.y();
    zValues(i) = qValue.z();
  }

  outValue.attr("class") = Rcpp::CharacterVector::create("qts", "tbl_df", "tbl", "data.frame");
  return outValue;
}

Rcpp::DataFrame derivative_qts_impl(const Rcpp::DataFrame &qts)
{
  unsigned int nGrid = qts.nrows();
  Rcpp::DataFrame outValue = Rcpp::clone(qts);
  Rcpp::NumericVector wValues = outValue["w"];
  Rcpp::NumericVector xValues = outValue["x"];
  Rcpp::NumericVector yValues = outValue["y"];
  Rcpp::NumericVector zValues = outValue["z"];
  Eigen::Quaterniond currentQValue, previousQvalue;
  currentQValue.w() = wValues(nGrid - 1);
  currentQValue.x() = xValues(nGrid - 1);
  currentQValue.y() = yValues(nGrid - 1);
  currentQValue.z() = zValues(nGrid - 1);

  for (unsigned int i = nGrid - 1;i > 0;--i)
  {
    previousQvalue.w() = wValues(i - 1);
    previousQvalue.x() = xValues(i - 1);
    previousQvalue.y() = yValues(i - 1);
    previousQvalue.z() = zValues(i - 1);

    currentQValue = previousQvalue.inverse() * currentQValue;

    wValues(i) = currentQValue.w();
    xValues(i) = currentQValue.x();
    yValues(i) = currentQValue.y();
    zValues(i) = currentQValue.z();

    currentQValue = previousQvalue;
  }

  return outValue;
}

Rcpp::DataFrame log_qts_impl(const Rcpp::DataFrame &qts)
{
  unsigned int nGrid = qts.nrows();
  Rcpp::DataFrame outValue = Rcpp::clone(qts);
  Rcpp::NumericVector wValues = outValue["w"];
  Rcpp::NumericVector xValues = outValue["x"];
  Rcpp::NumericVector yValues = outValue["y"];
  Rcpp::NumericVector zValues = outValue["z"];
  Eigen::Quaterniond qValue;

  for (unsigned int i = 0;i < nGrid;++i)
  {
    qValue = Eigen::Quaterniond(wValues(i), xValues(i), yValues(i), zValues(i));
    qValue = logq(qValue);
    wValues(i) = qValue.w();
    xValues(i) = qValue.x();
    yValues(i) = qValue.y();
    zValues(i) = qValue.z();
  }

  outValue.attr("class") = Rcpp::CharacterVector::create("qts", "tbl_df", "tbl", "data.frame");
  return outValue;
}

Rcpp::DataFrame exp_qts_impl(const Rcpp::DataFrame &qts)
{
  unsigned int nGrid = qts.nrows();
  Rcpp::DataFrame outValue = Rcpp::clone(qts);
  Rcpp::NumericVector wValues = outValue["w"];
  Rcpp::NumericVector xValues = outValue["x"];
  Rcpp::NumericVector yValues = outValue["y"];
  Rcpp::NumericVector zValues = outValue["z"];
  Eigen::Quaterniond qValue;

  for (unsigned int i = 0;i < nGrid;++i)
  {
    qValue = Eigen::Quaterniond(wValues(i), xValues(i), yValues(i), zValues(i));
    qValue = expq(qValue);
    wValues(i) = qValue.w();
    xValues(i) = qValue.x();
    yValues(i) = qValue.y();
    zValues(i) = qValue.z();
  }

  outValue.attr("class") = Rcpp::CharacterVector::create("qts", "tbl_df", "tbl", "data.frame");
  return outValue;
}

Rcpp::List centring_qts_impl(const Rcpp::DataFrame &qts, const bool standardize)
{
  unsigned int nGrid = qts.nrows();
  Rcpp::DataFrame outValue = Rcpp::clone(qts);
  Rcpp::NumericVector wValues = outValue["w"];
  Rcpp::NumericVector xValues = outValue["x"];
  Rcpp::NumericVector yValues = outValue["y"];
  Rcpp::NumericVector zValues = outValue["z"];

  std::vector<Eigen::VectorXd> qValues(nGrid);
  Eigen::Vector4d meanValue;
  for (unsigned int i = 0;i < nGrid;++i)
  {
    meanValue(0) = wValues(i);
    meanValue(1) = xValues(i);
    meanValue(2) = yValues(i);
    meanValue(3) = zValues(i);
    qValues[i] = meanValue;
  }

  meanValue = gmean(qValues);
  Eigen::Quaterniond meanQValue(meanValue(0), meanValue(1), meanValue(2), meanValue(3)), workQValue;
  meanQValue = meanQValue.inverse();

  for (unsigned int i = 0;i < nGrid;++i)
  {
    workQValue = Eigen::Quaterniond(wValues(i), xValues(i), yValues(i), zValues(i));
    workQValue = meanQValue * workQValue;
    wValues(i) = workQValue.w();
    xValues(i) = workQValue.x();
    yValues(i) = workQValue.y();
    zValues(i) = workQValue.z();
  }

  double sdValue = 0;
  if (standardize)
  {
    outValue = log_qts_impl(outValue);
    sdValue = std::sqrt(gvariance(qValues, meanValue));
    wValues = outValue["w"];
    xValues = outValue["x"];
    yValues = outValue["y"];
    zValues = outValue["z"];
    wValues = wValues / sdValue;
    xValues = xValues / sdValue;
    yValues = yValues / sdValue;
    zValues = zValues / sdValue;
    outValue = exp_qts_impl(outValue);
  }

  outValue.attr("class") = Rcpp::CharacterVector::create("qts", "tbl_df", "tbl", "data.frame");

  return Rcpp::List::create(
    Rcpp::Named("qts") = outValue,
    Rcpp::Named("mean") = meanValue,
    Rcpp::Named("sd") = sdValue
  );
}
