#ifndef EIGEN_H
#define EIGEN_H

#include <gsl/gsl_matrix.h>
#include "differential.h"
#include <utility>

namespace eigen
{
    double computeSpectralRadius();
    void init(const size_t size, const Differential::Type type, const double points[], const size_t npts);
    void freemem();
}

#endif // EIGEN_H
