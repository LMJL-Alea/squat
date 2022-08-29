// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_squat_RCPPEXPORTS_H_GEN_
#define RCPP_squat_RCPPEXPORTS_H_GEN_

#include "squat_types.h"
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

namespace squat {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("squat", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("squat", "_squat_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in squat");
            }
        }
    }

    inline Rcpp::DataFrame resample_qts(const Rcpp::DataFrame& qts, double tmin = NA_REAL, double tmax = NA_REAL, const unsigned int nout = 0, const bool disable_normalization = false) {
        typedef SEXP(*Ptr_resample_qts)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_resample_qts p_resample_qts = NULL;
        if (p_resample_qts == NULL) {
            validateSignature("Rcpp::DataFrame(*resample_qts)(const Rcpp::DataFrame&,double,double,const unsigned int,const bool)");
            p_resample_qts = (Ptr_resample_qts)R_GetCCallable("squat", "_squat_resample_qts");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_resample_qts(Shield<SEXP>(Rcpp::wrap(qts)), Shield<SEXP>(Rcpp::wrap(tmin)), Shield<SEXP>(Rcpp::wrap(tmax)), Shield<SEXP>(Rcpp::wrap(nout)), Shield<SEXP>(Rcpp::wrap(disable_normalization)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::DataFrame >(rcpp_result_gen);
    }

    inline Rcpp::DataFrame smooth_qts(const Rcpp::DataFrame& qts, const double alpha = 0.5) {
        typedef SEXP(*Ptr_smooth_qts)(SEXP,SEXP);
        static Ptr_smooth_qts p_smooth_qts = NULL;
        if (p_smooth_qts == NULL) {
            validateSignature("Rcpp::DataFrame(*smooth_qts)(const Rcpp::DataFrame&,const double)");
            p_smooth_qts = (Ptr_smooth_qts)R_GetCCallable("squat", "_squat_smooth_qts");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_smooth_qts(Shield<SEXP>(Rcpp::wrap(qts)), Shield<SEXP>(Rcpp::wrap(alpha)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::DataFrame >(rcpp_result_gen);
    }

    inline double GeodesicQuaternionDistance(const Rcpp::NumericMatrix& M1, const Rcpp::NumericMatrix& M2, const unsigned int index1, const unsigned int index2) {
        typedef SEXP(*Ptr_GeodesicQuaternionDistance)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_GeodesicQuaternionDistance p_GeodesicQuaternionDistance = NULL;
        if (p_GeodesicQuaternionDistance == NULL) {
            validateSignature("double(*GeodesicQuaternionDistance)(const Rcpp::NumericMatrix&,const Rcpp::NumericMatrix&,const unsigned int,const unsigned int)");
            p_GeodesicQuaternionDistance = (Ptr_GeodesicQuaternionDistance)R_GetCCallable("squat", "_squat_GeodesicQuaternionDistance");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_GeodesicQuaternionDistance(Shield<SEXP>(Rcpp::wrap(M1)), Shield<SEXP>(Rcpp::wrap(M2)), Shield<SEXP>(Rcpp::wrap(index1)), Shield<SEXP>(Rcpp::wrap(index2)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline Rcpp::NumericMatrix RegularizeGrid(const Rcpp::NumericVector& grid, const Rcpp::NumericMatrix& values, const double gridLowerBound, const double gridUpperBound, const unsigned int numberOfPoints) {
        typedef SEXP(*Ptr_RegularizeGrid)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_RegularizeGrid p_RegularizeGrid = NULL;
        if (p_RegularizeGrid == NULL) {
            validateSignature("Rcpp::NumericMatrix(*RegularizeGrid)(const Rcpp::NumericVector&,const Rcpp::NumericMatrix&,const double,const double,const unsigned int)");
            p_RegularizeGrid = (Ptr_RegularizeGrid)R_GetCCallable("squat", "_squat_RegularizeGrid");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_RegularizeGrid(Shield<SEXP>(Rcpp::wrap(grid)), Shield<SEXP>(Rcpp::wrap(values)), Shield<SEXP>(Rcpp::wrap(gridLowerBound)), Shield<SEXP>(Rcpp::wrap(gridUpperBound)), Shield<SEXP>(Rcpp::wrap(numberOfPoints)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::NumericMatrix >(rcpp_result_gen);
    }

    inline Rcpp::NumericMatrix GetGeodesicMean(const Rcpp::NumericMatrix& values) {
        typedef SEXP(*Ptr_GetGeodesicMean)(SEXP);
        static Ptr_GetGeodesicMean p_GetGeodesicMean = NULL;
        if (p_GetGeodesicMean == NULL) {
            validateSignature("Rcpp::NumericMatrix(*GetGeodesicMean)(const Rcpp::NumericMatrix&)");
            p_GetGeodesicMean = (Ptr_GetGeodesicMean)R_GetCCallable("squat", "_squat_GetGeodesicMean");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_GetGeodesicMean(Shield<SEXP>(Rcpp::wrap(values)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::NumericMatrix >(rcpp_result_gen);
    }

    inline Eigen::VectorXd gmean(const std::vector<Eigen::VectorXd>& quaternionSample, unsigned int maxIterations = 2000, double maxEpsilon = 1.0e-5) {
        typedef SEXP(*Ptr_gmean)(SEXP,SEXP,SEXP);
        static Ptr_gmean p_gmean = NULL;
        if (p_gmean == NULL) {
            validateSignature("Eigen::VectorXd(*gmean)(const std::vector<Eigen::VectorXd>&,unsigned int,double)");
            p_gmean = (Ptr_gmean)R_GetCCallable("squat", "_squat_gmean");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_gmean(Shield<SEXP>(Rcpp::wrap(quaternionSample)), Shield<SEXP>(Rcpp::wrap(maxIterations)), Shield<SEXP>(Rcpp::wrap(maxEpsilon)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Eigen::VectorXd >(rcpp_result_gen);
    }

    inline Eigen::VectorXd gmedian(const std::vector<Eigen::VectorXd>& quaternionSample, unsigned int maxIterations = 2000, double maxEpsilon = 1.0e-5) {
        typedef SEXP(*Ptr_gmedian)(SEXP,SEXP,SEXP);
        static Ptr_gmedian p_gmedian = NULL;
        if (p_gmedian == NULL) {
            validateSignature("Eigen::VectorXd(*gmedian)(const std::vector<Eigen::VectorXd>&,unsigned int,double)");
            p_gmedian = (Ptr_gmedian)R_GetCCallable("squat", "_squat_gmedian");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_gmedian(Shield<SEXP>(Rcpp::wrap(quaternionSample)), Shield<SEXP>(Rcpp::wrap(maxIterations)), Shield<SEXP>(Rcpp::wrap(maxEpsilon)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Eigen::VectorXd >(rcpp_result_gen);
    }

}

#endif // RCPP_squat_RCPPEXPORTS_H_GEN_
