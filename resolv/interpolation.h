#ifndef INTERPOLATION_H
#define INTERPOLATION_H


namespace interpolation
{
    void splineInit(const int size, const int sizem, const double mx[], const double my[], const double sx[], const double sy[]);
    double mu(const double x);
    double gamma(const double x);
    double beta(const double x1, const double x2);
    double S0(const double x);
    void free();
}

#endif // INTERPOLATION_H
