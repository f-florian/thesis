#ifndef EIGEN_H
#define EIGEN_H

#include <gsl/gsl_matrix.h>

namespace eigen
{
    double computeSpectralRadius();
    void init(const unsigned int size, double start, double end);
    void freemem();
}

#endif // EIGEN_H
