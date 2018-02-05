#include "differential.h"
#include <cmath>

const double pi=3.14159265358979323846264338327950288419716939937510;

Differential::Differential(unsigned short order_p)
{
	order=order_p;
	nodesx=new double(order);
	qw=new double(order);
	dw=new double(order);
	for(int i=0;i<order;i++){
	  nodesx[i]=cos(pi*(2*i+1)/(2*order));
	  //
	}
 	//calculate
}
double Differential::nodes(unsigned short index)
{
	//throw exception
	return nodesx[index];
}
double Differential::quadratureWeights(unsigned short index)
{
	return qw[index];
}
double Differential::differentiationWeights(unsigned short index)
{
	return dw[index];
}
