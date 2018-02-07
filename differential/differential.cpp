#include "differential.h"
#include <cmath>
#include <stdexcept>

const double pi=3.14159265358979323846264338327950288419716939937510;

Differential::Differential(unsigned short order_p)
{
  order=order_p;
	nodesx=new double(order);
	qw=new double(order);
	dw=new double(order*order);
	double difftmp[order*order];
	double prods[order];
	for(int i=0;i<order;i++){
	  nodesx[i]=cos(pi*(2*i+1)/(2*order));
	  prods[i]=1;
	  dw[i*order+i]=0;
	  // store differences in the lower triagular part
	  for(int j=0;j<i;j++){
	    difftmp[i*order+j]=nodesx[i]-nodesx[j];
	    prods[i]*=difftmp[i*order+j];
	    prods[j]*=-difftmp[i*order+j];
	    dw[i*order+i]+=1/difftmp[i*order+j];
	    dw[j*order+j]-=1/difftmp[i*order+j];
	  }
	}
	  for(int i=0;i<order;i++)
	    for(int j=0;j<i;j++) {
	      dw[i*order+j]=prods[j]/prods[i]/difftmp[i*order+j];
	      dw[j*order+i]=-prods[i]/prods[j]/difftmp[i*order+j];
	    }		
	  //check precomputed derivatives (manual)
	  //compute integrals using gsl Fixed point quadratures (p 185/199)
}
double Differential::nodes(unsigned short index)
{
  if(index>=order)
    throw(std::range_error("requested point not in matrix"));
  return nodesx[index];
}
double Differential::quadratureWeights(unsigned short index)
{
  if(index>=order)
    throw(std::range_error("requested point not in matrix"));
  return dw[index];
}
double Differential::differentiationWeights(unsigned short index, unsigned short point)
{
  if((index>=order)||(point>=order))
    throw(std::range_error("requested point not in matrix"));
  return qw[point*order+index];
}
