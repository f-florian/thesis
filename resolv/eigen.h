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
<<<<<<< HEAD
    void init(const size_t size, const Differential::Type type, const std::array<double,2> points);
=======
>>>>>>> 3a6e8bdefaf281450206d5f662e03107eb09711d
    void freemem();
}

#endif // EIGEN_H
