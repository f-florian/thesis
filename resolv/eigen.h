#ifndef EIGEN_H
#define EIGEN_H

#include <gsl/gsl_matrix.h>
#include "differential.h"
#include <utility>
#include <array>

namespace eigen
{
    double computeSpectralRadius(const double hint);
    void init(const size_t size, const Differential::Type type, const double points[], const size_t npts);
    template<size_t n> auto init(const size_t size, const Differential::Type type, const std::array<double,n> points) {return init(size,type,points.data(),n);}
    void freemem();
}

#endif // EIGEN_H
