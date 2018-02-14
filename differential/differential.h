#ifndef DIFFERENTIAL_H
#define DIFFERENTIAL_H

/**
 * store information on weights for interpolatory differential opertions (i.e. integral and derivative) with a given (minimum) order
 *
 * stores a mesh of [0,1] and weights wich can be used to compute \int_0^1 f(x) \de x and f'(y) for y a point of the mesh.
 * the (minimum) converdence order of the interpolating formulas obtained is specified when the class is constructed
 */
class Differential
{
public:
	Differential(unsigned short order);												//!<initialize data for a given order
	double nodes(unsigned short index);												//!< get index-th node in the mesh
	double quadratureWeights(unsigned short index);									//!< get index-th quadrature weight
	double differentiationWeights(unsigned short index, unsigned short point);		//!< get index-th weight for approximating derivative in point-th point
private:
	unsigned short order;
	double *nodesx;
	double *qw;
	double *dw;
};

#endif // DIFFERENTIAL_H
