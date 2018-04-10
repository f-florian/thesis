#ifndef EIGEN_H
#define EIGEN_H

#include <gsl/gsl_matrix.h>
#include "differential.h"
#include <utility>

namespace eigen
{
    double computeSpectralRadius();
    void init(const unsigned int size, double start, double end, Differential::Type type);
    void freemem();
}

#endif // EIGEN_H
