#include "functions.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <cmath>

namespace parameters{
    namespace {
        double length_m;
        constexpr double c0=1, D0=1, beta0=1, mu0=1;
        double c1, D1;
    }
    void init(const double length_p, const double c1_p, const double D1_p)
    {
        length_m=length_p;
        c1=c1_p;
        D1=D1_p;
    }
    double c(const double a)
    {
        return c0+c1*a;
    }
    double D(const double a)
    {
        return D0+D1*a;
    }
    double beta(const double a)
    {
        return beta0;
    }
    double mu(const double a);
    {
        return mu0;
    }
    double length()
    {
        return length_m;
    }
}
