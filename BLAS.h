#ifndef BLAS_H
#define BLAS_H

#include <complex>

#define srot srot_
#define drot drot_
#define crot crot_
#define zrot zrot_

#define BLAS_INT int

extern "C"
{
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
}

#endif
