#include "interpolation.h"

#include <RcppEigen.h>
#include <algorithm>

Rcpp::DataFrame resample_qts(const Rcpp::DataFrame &qts,
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

  double xmin = inputTimeValues(0);
  double xmax = inputTimeValues(sizeIn - 1);
  double step = (xmax - xmin) / (sizeOut - 1.0);

  auto posInf = inputTimeValues.begin();
  auto posSup = inputTimeValues.begin();

  Rcpp::NumericVector outputTimeValues = Rcpp::NumericVector(sizeOut);
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
    double newx = xmin + (double)i * step;
    outputTimeValues(i) = newx;

    if (*posInf == newx)
    {
      unsigned int idxInf = posInf - inputTimeValues.begin();
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
    else
    {
      auto oldPosInf = posInf;
      auto oldPosSup = posSup;

      if (*posSup <= newx)
      {
        posSup = std::upper_bound(posInf, inputTimeValues.end(), newx);
        if (posSup == inputTimeValues.end())
        {
          posSup = inputTimeValues.end() - 1;
          posInf = inputTimeValues.end() - 1;
        }
        else
          posInf = posSup - 1;
      }

      unsigned int idxInf = posInf - inputTimeValues.begin();
      unsigned int idxSup = posSup - inputTimeValues.begin();

      if (posInf != oldPosInf)
      {
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

      if (posSup != oldPosSup)
        xsup = *posSup;

      if (xsup > xinf)
      {
        double range = xsup - xinf;
        double alpha = (newx - xinf) / range;

        if (posSup != oldPosSup)
        {
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

        Qinf = Qinf.slerp(alpha, Qsup);
      }
    }

    outputWValues(i) = Qinf.w();
    outputXValues(i) = Qinf.x();
    outputYValues(i) = Qinf.y();
    outputZValues(i) = Qinf.z();
  }

  return outputValue;
}

