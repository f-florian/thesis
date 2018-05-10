#include "functions.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <cmath>

namespace parameters{
    namespace {
        double length_m;
        double c0, D0, beta0, mu0, c1, D1, p;
    }
    void init(const double length_p, const double c0_p, const double D0_p, const double beta0_p, const double mu0_p, const double c1_p, const double D1_p, const int p_p)
    {
        length_m=length_p;
        c1=c1_p;
        D1=D1_p;
        c0=c0_p;
        D0=D0_p;
        beta0=beta0_p;
        mu0=mu0_p;
        p=p_p;
    }
    double c(const double a)
    {
        double s=c0+c1*pow(a,p);
        // printf("c(%e) %e\n", a, s);
        return s;
    }
    double D(const double a)
    {
        double s=D0+D1*pow(a,p);
        // printf("D(%e) %e\n", a, s);
        return s;
    }
    double beta([[maybe_unused]] const double a)
    {
        // return beta0*a*(length_m-a);
        return beta0;
    }
    double betaonc(const double a)
    {
        return beta0*a/c0;
    }
    double mu([[maybe_unused]] const double a)
    {
        return mu0;
    }
    double length()
    {
        return length_m;
    }
}
