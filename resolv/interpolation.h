#ifndef INTERPOLATION_H
#define INTERPOLATION_H


namespace interpolation
{
    void splineInit(int size, int sizem, const double mx[], const double my[], const double sx[], const double sy[]);
    double mu(double x);
    double gamma(double x);
    double beta(double x1, double x2);
    double S0(double x);
}

#endif // INTERPOLATION_H
