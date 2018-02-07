#ifndef DIFFERENTIAL_H
#define DIFFERENTIAL_H


class Differential
{
public:
	Differential(unsigned short order);
	double nodes(unsigned short index);
	double quadratureWeights(unsigned short index);
	double differentiationWeights(unsigned short index, unsigned short point);
private:
	unsigned short order;
	double *nodesx;
	double *qw;
	double *dw;
};

#endif // DIFFERENTIAL_H
