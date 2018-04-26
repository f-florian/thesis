#ifndef EIGEN_H
#define EIGEN_H

#include <gsl/gsl_matrix.h>
#include "differential.h"
#include <utility>
#include <array>

namespace eigen
{
    double computeSpectralRadius();
    void init(const size_t size, const Differential::Type type, const double points[], const size_t npts);
    void init(const size_t size, const Differential::Type type, const std::array<double,2> points);
    void freemem();
}

#endif // EIGEN_H
