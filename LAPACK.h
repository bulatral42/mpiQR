#ifndef LAPACK_H
#define LAPACK_H

#define slartg slartg_
#define dlartg dlartg_
#define clartg clartg_
#define zlartg zlartg_

#define LAPACK_INT int

extern "C"
{
    void slartg(const float &, const float &, float &, float &, float &);

    void dlartg(const double &, const double &, double &, double &, double &);

    void clartg(const std::complex<float> &, const std::complex<float> &, float &, std::complex<float> &, std::complex<float> &);

    void zlartg(const std::complex<double> &, const std::complex<double> &, double &, std::complex<double> &, std::complex<double> &);
}

#endif
