#ifndef DIFFERENTIAL_H
#define DIFFERENTIAL_H

/**
 * store information on weights for interpolatory differential opertions (i.e. integral and derivative) with a given (minimum) order
 *
 * stores a mesh of [0,1] and weights wich can be used to compute \int_0^1 f(x) \de x and f'(y) for y a point of the mesh.
 * the (minimum) convergence order of the interpolating formulas obtained is specified when the class is constructed
 */
class Differential
{
public:
    Differential(unsigned short order);                                                                             //!<initialize data for a given order
    double nodes(unsigned short index, double end=1, double start=0);                                               //!< get index-th node in the mesh for the interval [start-end]
    double quadratureWeights(unsigned short index, double end=1, double start=0);	                                //!< get index-th quadrature weight, properly scaled for nodes in [start, end] 
    double differentiationWeights(unsigned short index, unsigned short point, double end=1, double start=0);        //!< get index-th weight for approximating derivative in point-th point, properly scaled for nodes in [start, end]
    double* StealNodes();
    double* StealQuadratureWeights();
    double* StealDifferentiationWeights();
private:
    unsigned short order;
    double *nodesx;
    double *qw;
    double *dw;
};

#endif // DIFFERENTIAL_H
