#ifndef BLAS_H
#define BLAS_H

#include <complex>

#define sswap sswap_
#define dswap dswap_
#define cswap cswap_
#define zswap zswap_
#define saxpy saxpy_
#define daxpy daxpy_
#define caxpy caxpy_
#define zaxpy zaxpy_
#define scopy scopy_
#define dcopy dcopy_
#define ccopy ccopy_
#define zcopy zcopy_
#define snrm2 snrm2_
#define dnrm2 dnrm2_
#define scnrm2 scnrm2_
#define dznrm2 dznrm2_
#define sdot sdot_
#define ddot ddot_
#define cdotc_stackreturn cdotc_stackreturn_
#define cdotu_stackreturn cdotu_stackreturn_
#define zdotc_stackreturn zdotc_stackreturn_
#define zdotu_stackreturn zdotu_stackreturn_
#define sscal sscal_
#define dscal dscal_
#define cscal cscal_
#define csscal csscal_
#define zscal zscal_
#define zdscal zdscal_
#define srot srot_
#define drot drot_
#define crot crot_
#define zrot zrot_
#define sger sger_
#define dger dger_
#define cgerc cgerc_
#define zgerc zgerc_
#define cgeru cgeru_
#define zgeru zgeru_
#define sgemv sgemv_
#define dgemv dgemv_
#define cgemv cgemv_
#define zgemv zgemv_
#define stpsv stpsv_
#define dtpsv dtpsv_
#define ctpsv ctpsv_
#define ztpsv ztpsv_
#define sgemm sgemm_
#define dgemm dgemm_
#define cgemm cgemm_
#define zgemm zgemm_

#define BLAS_INT int

extern "C"
{
    void sswap(const BLAS_INT &,
               const float *, const BLAS_INT &,
               float *, const BLAS_INT &);

    void dswap(const BLAS_INT &,
               const double *, const BLAS_INT &,
               double *, const BLAS_INT &);

    void cswap(const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &,
               std::complex<float> *, const BLAS_INT &);

    void zswap(const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &,
               std::complex<double> *, const BLAS_INT &);

    void saxpy(const BLAS_INT &,
               const float &,
               const float *, const BLAS_INT &,
               float *, const BLAS_INT &);

    void daxpy(const BLAS_INT &,
               const double &,
               const double *, const BLAS_INT &,
               double *, const BLAS_INT &);

    void caxpy(const BLAS_INT &,
               const std::complex<float> &,
               const std::complex<float> *, const BLAS_INT &,
               std::complex<float> *, const BLAS_INT &);

    void zaxpy(const BLAS_INT &,
               const std::complex<double> &,
               const std::complex<double> *, const BLAS_INT &,
               std::complex<double> *, const BLAS_INT &);

    void scopy(const BLAS_INT &,
               const float *, const BLAS_INT &,
               float *, const BLAS_INT &);

    void dcopy(const BLAS_INT &,
               const double *, const BLAS_INT &,
               double *, const BLAS_INT &);

    void ccopy(const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &,
               std::complex<float> *, const BLAS_INT &);

    void zcopy(const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &,
               std::complex<double> *, const BLAS_INT &);

    float sdot(const BLAS_INT &,
               const float *, const BLAS_INT &,
               const float *, const BLAS_INT &);

    double ddot(const BLAS_INT &,
               const double *, const BLAS_INT &,
               const double *, const BLAS_INT &);

    void cdotu_stackreturn(std::complex<float> &, const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &);

    void cdotc_stackreturn(std::complex<float> &, const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &);

    void zdotu_stackreturn(std::complex<double> &, const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &);

    void zdotc_stackreturn(std::complex<double> &, const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &);

    float snrm2(const BLAS_INT &,
                const float *, const BLAS_INT &);

    double dnrm2(const BLAS_INT &,
                  const double *, const BLAS_INT &);

    float scnrm2(const BLAS_INT &,
                 const std::complex<float> *, const BLAS_INT &);

    double dznrm2(const BLAS_INT &,
                  const std::complex<double> *, const BLAS_INT &);

    void sscal(const BLAS_INT &,
               const float &,
               const float *, const BLAS_INT &);

    void dscal(const BLAS_INT &,
               const double &,
               const double *, const BLAS_INT &);

    void cscal(const BLAS_INT &,
               const std::complex<float> &,
               const std::complex<float> *, const BLAS_INT &);

    void csscal(const BLAS_INT &,
                const float &,
                const std::complex<float> *, const BLAS_INT &);

    void zscal(const BLAS_INT &,
               const std::complex<double> &,
               const std::complex<double> *, const BLAS_INT &);

    void zdscal(const BLAS_INT &,
                const double &,
                const std::complex<double> *, const BLAS_INT &);

    void srot(const BLAS_INT &,
              float *, const BLAS_INT &,
              float *, const BLAS_INT &,
              const float &, const float &);

    void drot(const BLAS_INT &,
              double *, const BLAS_INT &,
              double *, const BLAS_INT &,
              const double &, const double &);

    void crot(const BLAS_INT &,
              std::complex<float> *, const BLAS_INT &,
              std::complex<float> *, const BLAS_INT &,
              const float &, const std::complex<float> &);

    void zrot(const BLAS_INT &,
              std::complex<double> *, const BLAS_INT &,
              std::complex<double> *, const BLAS_INT &,
              const double &, const std::complex<double> &);

    void sger(const BLAS_INT &, const BLAS_INT &,
              const float &,
              const float *, const BLAS_INT &,
              const float *, const BLAS_INT &,
              float *, const BLAS_INT &);

    void dger(const BLAS_INT &, const BLAS_INT &,
              const double &,
              const double *, const BLAS_INT &,
              const double *, const BLAS_INT &,
              double *, const BLAS_INT &);

    void cgerc(const BLAS_INT &, const BLAS_INT &,
               const std::complex<float> &,
               const std::complex<float> *, const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &,
               std::complex<float> *, const BLAS_INT &);

    void zgerc(const BLAS_INT &, const BLAS_INT &,
               const std::complex<double> &,
               const std::complex<double> *, const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &,
               std::complex<double> *, const BLAS_INT &);

    void cgeru(const BLAS_INT &, const BLAS_INT &,
               const std::complex<float> &,
               const std::complex<float> *, const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &,
               std::complex<float> *, const BLAS_INT &);

    void zgeru(const BLAS_INT &, const BLAS_INT &,
               const std::complex<double> &,
               const std::complex<double> *, const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &,
               std::complex<double> *, const BLAS_INT &);

    void sgemv(const char &,
               const BLAS_INT &, const BLAS_INT &,
               const float &,
               const float *, const BLAS_INT &,
               const float *, const BLAS_INT &,
               const float &, 
               float *, const BLAS_INT &);

    void dgemv(const char &,
               const BLAS_INT &, const BLAS_INT &,
               const double &,
               const double *, const BLAS_INT &,
               const double *, const BLAS_INT &,
               const double &, 
               double *, const BLAS_INT &);

    void cgemv(const char &,
               const BLAS_INT &, const BLAS_INT &,
               const std::complex<float> &,
               const std::complex<float> *, const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &,
               const std::complex<float> &, 
               std::complex<float> *, const BLAS_INT &);

    void zgemv(const char &,
               const BLAS_INT &, const BLAS_INT &,
               const std::complex<double> &,
               const std::complex<double> *, const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &,
               const std::complex<double> &, 
               std::complex<double> *, const BLAS_INT &);

    void stpsv(const char &,
               const char &,
               const char &,
               const BLAS_INT &,
               const float *,
               float *, const BLAS_INT &);

    void dtpsv(const char &,
               const char &,
               const char &,
               const BLAS_INT &,
               const double *,
               double *, const BLAS_INT &);

    void ctpsv(const char &,
               const char &,
               const char &,
               const BLAS_INT &,
               const std::complex<float> *,
               std::complex<float> *, const BLAS_INT &);

    void ztpsv(const char &,
               const char &,
               const char &,
               const BLAS_INT &,
               const std::complex<double> *,
               std::complex<double> *, const BLAS_INT &);

    void sgemm(const char &,
               const char &,
               const BLAS_INT &, const BLAS_INT &,
               const BLAS_INT &,
               const float &,
               const float *, const BLAS_INT &,
               const float *, const BLAS_INT &,
               const float &,
               float *, const BLAS_INT &);

    void dgemm(const char &,
               const char &,
               const BLAS_INT &, const BLAS_INT &,
               const BLAS_INT &,
               const double &,
               const double *, const BLAS_INT &,
               const double *, const BLAS_INT &,
               const double &,
               double *, const BLAS_INT &);

    void cgemm(const char &,
               const char &,
               const BLAS_INT &, const BLAS_INT &,
               const BLAS_INT &,
               const std::complex<float> &,
               const std::complex<float> *, const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &,
               const std::complex<float> &,
               std::complex<float> *, const BLAS_INT &);

    void zgemm(const char &,
               const char &,
               const BLAS_INT &, const BLAS_INT &,
               const BLAS_INT &,
               const std::complex<double> &,
               const std::complex<double> *, const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &,
               const std::complex<double> &,
               std::complex<double> *, const BLAS_INT &);
}

#endif
