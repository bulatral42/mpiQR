#ifndef LAPACK_H
#define LAPACK_H

#define slartg slartg_
#define dlartg dlartg_
#define clartg clartg_
#define zlartg zlartg_
#define clacgv clacgv_
#define zlacgv zlacgv_
#define stptrs stptrs_
#define dtptrs dtptrs_
#define ctptrs ctptrs_
#define ztptrs ztptrs_

#define LAPACK_INT int

extern "C"
{
    void slartg(const float &, const float &, float &, float &, float &);

    void dlartg(const double &, const double &, double &, double &, double &);

    void clartg(const std::complex<float> &, const std::complex<float> &, float &, std::complex<float> &, std::complex<float> &);

    void zlartg(const std::complex<double> &, const std::complex<double> &, double &, std::complex<double> &, std::complex<double> &);

    void clacgv(const LAPACK_INT &, std::complex<float> *, const LAPACK_INT &);

    void zlacgv(const LAPACK_INT &, std::complex<double> *, const LAPACK_INT &);

    void stptrs(const char &, const char &, const char &, const LAPACK_INT &, const LAPACK_INT &,
                const float *, const float *, const LAPACK_INT &, const LAPACK_INT &);

    void dtptrs(const char &, const char &, const char &, const LAPACK_INT &, const LAPACK_INT &,
                const double *, const double *, const LAPACK_INT &, const LAPACK_INT &);

    void ctptrs(const char &, const char &, const char &, const LAPACK_INT &, const LAPACK_INT &,
                const std::complex<float> *, const std::complex<float> *, const LAPACK_INT &, const LAPACK_INT &);

    void ztptrs(const char &, const char &, const char &, const LAPACK_INT &, const LAPACK_INT &,
                const std::complex<double> *, const std::complex<double> *, const LAPACK_INT &, const LAPACK_INT &);
}

#endif
