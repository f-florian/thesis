#ifndef EIGEN_H
#define EIGEN_H

#include <gsl/gsl_matrix.h>

namespace eigen
{
    double computeSpectralRadius();
	void init(const unsigned int size, double start, double delta, unsigned short order);
    void freemem();
	void setOrder(unsigned short ni);
}

#endif // EIGEN_H
