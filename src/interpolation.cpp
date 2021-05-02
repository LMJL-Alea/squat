#include "interpolation.h"

#include <RcppEigen.h>
#include <algorithm>

Rcpp::DataFrame resample_qts(const Rcpp::DataFrame &qts,
                             double tmin,
                             double tmax,
                             const unsigned int nout,
                             const bool disable_normalization)
{
  // Assumes qts$time is sorted in ascending order
  int sizeIn = qts.nrows();
  unsigned int sizeOut = (nout == 0) ? sizeIn : nout;

  Eigen::Quaterniond Qinf, Qsup;
  double xinf, xsup;

  Rcpp::NumericVector inputTimeValues = qts["time"];
  Rcpp::NumericVector inputWValues = qts["w"];
  Rcpp::NumericVector inputXValues = qts["x"];
  Rcpp::NumericVector inputYValues = qts["y"];
  Rcpp::NumericVector inputZValues = qts["z"];

  if (R_IsNA(tmin))
    tmin = inputTimeValues(0);
  if (R_IsNA(tmax))
    tmax = inputTimeValues(sizeIn - 1);

  auto posInf = inputTimeValues.begin();
  auto posSup = inputTimeValues.begin();
  auto oldPosInf = inputTimeValues.end();
  auto oldPosSup = inputTimeValues.end();

  Rcpp::NumericVector outputTimeValues = Rcpp::wrap(Eigen::VectorXd::LinSpaced(sizeOut, tmin, tmax));;
  Rcpp::NumericVector outputWValues = Rcpp::NumericVector(sizeOut);
  Rcpp::NumericVector outputXValues = Rcpp::NumericVector(sizeOut);
  Rcpp::NumericVector outputYValues = Rcpp::NumericVector(sizeOut);
  Rcpp::NumericVector outputZValues = Rcpp::NumericVector(sizeOut);

  Rcpp::DataFrame outputValue = Rcpp::DataFrame::create(
    Rcpp::Named("time") = outputTimeValues,
    Rcpp::Named("w") = outputWValues,
    Rcpp::Named("x") = outputXValues,
    Rcpp::Named("y") = outputYValues,
    Rcpp::Named("z") = outputZValues
  );
  outputValue.attr("class") = Rcpp::CharacterVector::create("tbl_df", "tbl", "data.frame");

  double epsValue = std::sqrt(std::numeric_limits<double>::epsilon());
  Rcpp::LogicalVector isNormalized = Rcpp::rep(false, sizeIn);

  for (unsigned int i = 0;i < sizeOut;++i)
  {
    // Grab new time point
    double tnew = outputTimeValues(i);

    // Assign NA quaternion to output if new point is not in QTS range
    if (tnew < inputTimeValues(0) || tnew > inputTimeValues(sizeIn - 1))
    {
      outputWValues(i) = NA_REAL;
      outputXValues(i) = NA_REAL;
      outputYValues(i) = NA_REAL;
      outputZValues(i) = NA_REAL;
      continue;
    }

    if (posSup == inputTimeValues.end())
    {
      outputWValues(i) = inputTimeValues(sizeIn - 1);
      outputXValues(i) = inputTimeValues(sizeIn - 1);
      outputYValues(i) = inputTimeValues(sizeIn - 1);
      outputZValues(i) = inputTimeValues(sizeIn - 1);
      continue;
    }

    posSup = std::upper_bound(posSup, inputTimeValues.end(), tnew);
    posInf = posSup - 1;

    if (posInf != oldPosInf)
    {
      unsigned int idxInf = posInf - inputTimeValues.begin();
      xinf = *posInf;
      Qinf.w() = inputWValues(idxInf);
      Qinf.x() = inputXValues(idxInf);
      Qinf.y() = inputYValues(idxInf);
      Qinf.z() = inputZValues(idxInf);
      if (!isNormalized(idxInf) && !disable_normalization)
      {
        Qinf.normalize();
        inputWValues(idxInf) = Qinf.w();
        inputXValues(idxInf) = Qinf.x();
        inputYValues(idxInf) = Qinf.y();
        inputZValues(idxInf) = Qinf.z();
        isNormalized(idxInf) = true;
      }
    }

    if (posSup != inputTimeValues.end())
    {
      if (posSup != oldPosSup)
      {
        unsigned int idxSup = posSup - inputTimeValues.begin();
        xsup = *posSup;
        Qsup.w() = inputWValues(idxSup);
        Qsup.x() = inputXValues(idxSup);
        Qsup.y() = inputYValues(idxSup);
        Qsup.z() = inputZValues(idxSup);
        if (!isNormalized(idxSup) && !disable_normalization)
        {
          Qsup.normalize();
          inputWValues(idxSup) = Qsup.w();
          inputXValues(idxSup) = Qsup.x();
          inputYValues(idxSup) = Qsup.y();
          inputZValues(idxSup) = Qsup.z();
          isNormalized(idxSup) = true;
        }
      }

      if (xsup > xinf)
      {
        double range = xsup - xinf;
        double alpha = (tnew - xinf) / range;
        Qinf = Qinf.slerp(alpha, Qsup);
      }
    }

    outputWValues(i) = Qinf.w();
    outputXValues(i) = Qinf.x();
    outputYValues(i) = Qinf.y();
    outputZValues(i) = Qinf.z();

    oldPosSup = posSup;
    oldPosInf = posInf;
  }

  return outputValue;
}

Rcpp::DataFrame smooth_qts(const Rcpp::DataFrame &qts,
                           const double alpha)
{
  unsigned int nGrid = qts.nrows();
  Rcpp::DataFrame outValue = Rcpp::clone(qts);
  Rcpp::NumericVector wValues = outValue["w"];
  Rcpp::NumericVector xValues = outValue["x"];
  Rcpp::NumericVector yValues = outValue["y"];
  Rcpp::NumericVector zValues = outValue["z"];
  std::vector<Eigen::Quaterniond> qValues(nGrid);

  for (unsigned int i = 0;i < nGrid;++i)
  {
    qValues[i] = Eigen::Quaterniond(wValues(i), xValues(i), yValues(i), zValues(i));
    if (i == 0)
      continue;
    qValues[i] = qValues[i].slerp(alpha, qValues[i - 1]);
  }

  for (int i = nGrid - 2;i >= 0;--i)
  {
    qValues[i] = qValues[i].slerp(alpha, qValues[i + 1]);
    wValues(i) = qValues[i].w();
    xValues(i) = qValues[i].x();
    yValues(i) = qValues[i].y();
    zValues(i) = qValues[i].z();
  }

  return outValue;
}
