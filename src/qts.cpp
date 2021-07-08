#include "qts.h"
#include "rotations.h"
#include "representations.h"
#include <RcppEigen.h>

Rcpp::DataFrame qts2distance(const Rcpp::DataFrame &first_qts,
                             const Rcpp::DataFrame &second_qts)
{
  unsigned int nGrid = first_qts.nrows();
  Rcpp::NumericVector firstWValues = first_qts["w"];
  Rcpp::NumericVector firstXValues = first_qts["x"];
  Rcpp::NumericVector firstYValues = first_qts["y"];
  Rcpp::NumericVector firstZValues = first_qts["z"];
  Rcpp::NumericVector secondWValues = second_qts["w"];
  Rcpp::NumericVector secondXValues = second_qts["x"];
  Rcpp::NumericVector secondYValues = second_qts["y"];
  Rcpp::NumericVector secondZValues = second_qts["z"];
  Eigen::Quaterniond firstQValue, secondQValue;

  Rcpp::NumericVector distanceValues(nGrid);
  for (unsigned int i = 0;i < nGrid;++i)
  {
    firstQValue = Eigen::Quaterniond(
      firstWValues(i),
      firstXValues(i),
      firstYValues(i),
      firstZValues(i)
    );
    secondQValue = Eigen::Quaterniond(
      secondWValues(i),
      secondXValues(i),
      secondYValues(i),
      secondZValues(i)
    );
    distanceValues(i) = secondQValue.angularDistance(firstQValue);
  }

  Rcpp::DataFrame outValue = Rcpp::DataFrame::create(
    Rcpp::Named("time") = first_qts["time"],
    Rcpp::Named("distance") = distanceValues
  );

  outValue.attr("class") = Rcpp::CharacterVector::create("tbl_df", "tbl", "data.frame");

  return outValue;
}

Rcpp::DataFrame qts2norm(const Rcpp::DataFrame &qts,
                         const bool disable_normalization)
{
  unsigned int nSamples = qts.nrows();
  Eigen::Quaterniond qValue;
  Rcpp::NumericVector normValues(nSamples);
  Rcpp::NumericVector wValues = qts["w"];
  Rcpp::NumericVector xValues = qts["x"];
  Rcpp::NumericVector yValues = qts["y"];
  Rcpp::NumericVector zValues = qts["z"];

  Eigen::Quaterniond refValue;
  refValue.w() = 1.0;
  refValue.x() = 0.0;
  refValue.y() = 0.0;
  refValue.z() = 0.0;
  if (!disable_normalization)
    refValue.normalize();

  for (unsigned int i = 0;i < nSamples;++i)
  {
    qValue.w() = wValues(i);
    qValue.x() = xValues(i);
    qValue.y() = yValues(i);
    qValue.z() = zValues(i);
    if (!disable_normalization)
      qValue.normalize();

    normValues(i) = qValue.angularDistance(refValue);
  }

  Rcpp::DataFrame outValue = Rcpp::DataFrame::create(
    Rcpp::Named("time") = qts["time"],
    Rcpp::Named("norm") = normValues
  );

  outValue.attr("class") = Rcpp::CharacterVector::create("tbl_df", "tbl", "data.frame");

  return outValue;
}

Rcpp::DataFrame qts2angle(const Rcpp::DataFrame &qts,
                          const bool disable_normalization)
{
  unsigned int nSamples = qts.nrows();
  Eigen::Quaterniond qValue;
  Rcpp::NumericVector angleValues(nSamples);
  Rcpp::NumericVector wValues = qts["w"];
  Rcpp::NumericVector xValues = qts["x"];
  Rcpp::NumericVector yValues = qts["y"];
  Rcpp::NumericVector zValues = qts["z"];

  Eigen::Quaterniond refValue;
  refValue.w() = wValues(0);
  refValue.x() = xValues(0);
  refValue.y() = yValues(0);
  refValue.z() = zValues(0);
  if (!disable_normalization)
    refValue.normalize();

  for (unsigned int i = 0;i < nSamples;++i)
  {
    qValue.w() = wValues(i);
    qValue.x() = xValues(i);
    qValue.y() = yValues(i);
    qValue.z() = zValues(i);
    if (!disable_normalization)
      qValue.normalize();

    angleValues(i) = qValue.angularDistance(refValue);
  }

  Rcpp::DataFrame outValue = Rcpp::DataFrame::create(
    Rcpp::Named("time") = qts["time"],
    Rcpp::Named("angle") = angleValues
  );

  outValue.attr("class") = Rcpp::CharacterVector::create("tbl_df", "tbl", "data.frame");

  return outValue;
}

Rcpp::DataFrame reorient_qts(const Rcpp::DataFrame &qts,
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

  return resValue;
}

Rcpp::DataFrame normalize_qts(const Rcpp::DataFrame &qts)
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

Rcpp::DataFrame log_qts(const Rcpp::DataFrame &qts)
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

  return outValue;
}

Rcpp::DataFrame exp_qts(const Rcpp::DataFrame &qts)
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

  return outValue;
}

Rcpp::DataFrame centring_qts(const Rcpp::DataFrame &qts, const bool standardize)
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

  if (standardize)
  {
    outValue = log_qts(outValue);
    double sdValue = std::sqrt(gvariance(qValues, meanValue));
    wValues = outValue["w"];
    xValues = outValue["x"];
    yValues = outValue["y"];
    zValues = outValue["z"];
    wValues = wValues / sdValue;
    xValues = xValues / sdValue;
    yValues = yValues / sdValue;
    zValues = zValues / sdValue;
    outValue = exp_qts(outValue);
  }

  return outValue;
}

Rcpp::DataFrame mean_qts(const Rcpp::List &qts_list)
{
  unsigned int nSamples = qts_list.size();
  Rcpp::DataFrame outValue, tmpValue;
  outValue = Rcpp::clone(qts_list)[0];
  unsigned int nGrid = outValue.nrows();
  Rcpp::NumericVector wValues, xValues, yValues, zValues;
  std::vector<Eigen::VectorXd> qValues(nSamples);
  Eigen::Vector4d avgQValue;

  for (unsigned int i = 0;i < nGrid;++i)
  {
    for (unsigned int j = 0;j < nSamples;++j)
    {
      tmpValue = qts_list[j];
      wValues = tmpValue["w"];
      xValues = tmpValue["x"];
      yValues = tmpValue["y"];
      zValues = tmpValue["z"];
      avgQValue = {wValues(i), xValues(i), yValues(i), zValues(i)};
      qValues[j] = avgQValue;
    }

    avgQValue = gmean(qValues);

    wValues = outValue["w"];
    xValues = outValue["x"];
    yValues = outValue["y"];
    zValues = outValue["z"];
    wValues(i) = avgQValue(0);
    xValues(i) = avgQValue(1);
    yValues(i) = avgQValue(2);
    zValues(i) = avgQValue(3);
  }

  return outValue;
}

Rcpp::DataFrame median_qts(const Rcpp::List &qts_list)
{
  unsigned int nSamples = qts_list.size();
  Rcpp::DataFrame outValue, tmpValue;
  outValue = Rcpp::clone(qts_list)[0];
  unsigned int nGrid = outValue.nrows();
  Rcpp::NumericVector wValues, xValues, yValues, zValues;
  std::vector<Eigen::VectorXd> qValues(nSamples);
  Eigen::Vector4d avgQValue;

  for (unsigned int i = 0;i < nGrid;++i)
  {
    for (unsigned int j = 0;j < nSamples;++j)
    {
      tmpValue = qts_list[j];
      wValues = tmpValue["w"];
      xValues = tmpValue["x"];
      yValues = tmpValue["y"];
      zValues = tmpValue["z"];
      avgQValue = {wValues(i), xValues(i), yValues(i), zValues(i)};
      qValues[j] = avgQValue;
    }

    avgQValue = gmedian(qValues);

    wValues = outValue["w"];
    xValues = outValue["x"];
    yValues = outValue["y"];
    zValues = outValue["z"];
    wValues(i) = avgQValue(0);
    xValues(i) = avgQValue(1);
    yValues(i) = avgQValue(2);
    zValues(i) = avgQValue(3);
  }

  return outValue;
}

Rcpp::DataFrame qts2avts(const Rcpp::DataFrame &qts, const Rcpp::String &fixed_frame)
{
  unsigned int nGrid = qts.nrows();
  Rcpp::NumericVector inputTValues = qts["time"];
  Rcpp::NumericVector inputWValues = qts["w"];
  Rcpp::NumericVector inputXValues = qts["x"];
  Rcpp::NumericVector inputYValues = qts["y"];
  Rcpp::NumericVector inputZValues = qts["z"];
  Rcpp::NumericVector outputTValues(nGrid - 1);
  Rcpp::NumericVector outputXValues(nGrid - 1);
  Rcpp::NumericVector outputYValues(nGrid - 1);
  Rcpp::NumericVector outputZValues(nGrid - 1);

  Eigen::Quaterniond prevQValue, currQValue;
  for (unsigned int i = 1;i < nGrid;++i)
  {
    prevQValue = Eigen::Quaterniond(inputWValues(i - 1), inputXValues(i - 1), inputYValues(i - 1), inputZValues(i - 1));
    currQValue = Eigen::Quaterniond(inputWValues(i), inputXValues(i), inputYValues(i), inputZValues(i));
    double deltaTime = inputTValues(i) - inputTValues(i - 1);

    if (fixed_frame == "global")
      currQValue = currQValue * prevQValue.inverse();
    else if (fixed_frame == "body")
      currQValue = prevQValue.inverse() * currQValue;

    currQValue.coeffs() *= (2.0 / deltaTime);
    outputTValues(i - 1) = inputTValues(i);
    outputXValues(i - 1) = currQValue.x();
    outputYValues(i - 1) = currQValue.y();
    outputZValues(i - 1) = currQValue.z();
  }

  Rcpp::DataFrame outValue = Rcpp::DataFrame::create(
    Rcpp::Named("time") = outputTValues,
    Rcpp::Named("x") = outputXValues,
    Rcpp::Named("y") = outputYValues,
    Rcpp::Named("z") = outputZValues
  );

  outValue.attr("class") = Rcpp::CharacterVector::create("tbl_df", "tbl", "data.frame");

  return outValue;
}

Rcpp::DataFrame avts2qts(const Rcpp::DataFrame &avts,
                         const double init_t,
                         const Rcpp::NumericVector init_q)
{
  Rcpp::DataFrame outValue;
  return outValue;
}
