#include <cmath>
#include <stdexcept>
#include <gsl/gsl_integration.h>
#include "differential.h"

const double pi=3.14159265358979323846264338327950288419716939937510;

Differential::Differential(unsigned short order_p)
{
  order=order_p;
  dw=new double[order*order];
  gsl_integration_fixed_workspace *iws=gsl_integration_fixed_alloc(gsl_integration_fixed_legendre, order, 0,1,0,0);
  nodesx=iws->x;
  iws->x=NULL;
  qw=iws->weights;
  iws->weights=NULL;
  gsl_integration_fixed_free(iws);
  // for(int i=0;i<order;i++)
    // printf("%d %le %le\n",i, nodesx[i], qw[i]);
  double difftmp[order*order];
  double prods[order];
  for(int i=0;i<order;i++){
    prods[i]=1;
    dw[i*order+i]=0;
    // store differences in the lower triagular part; todo: only alloc lower triangular
    for(int j=0;j<i;j++){
      difftmp[i*order+j]=nodesx[i]-nodesx[j];
      prods[i]*=difftmp[i*order+j];
      prods[j]*=-difftmp[i*order+j];
      dw[i*order+i]+=1/difftmp[i*order+j];
      dw[j*order+j]-=1/difftmp[i*order+j];
    }
  }
  // for(int i=0;i<order*order;i++)
  // 	printf("dw[%d]=%lf\n",i, dw[i]);
  for(int i=0;i<order;i++)
    for(int j=0;j<i;j++) {
      // li'(xj)
      dw[i*order+j]=-prods[j]/prods[i]/difftmp[i*order+j];
      // lj'(xi)
      dw[j*order+i]=prods[i]/prods[j]/difftmp[i*order+j];
      // printf("set (%d,%d) (=%d,%d) to %lf (and %lf)\n",i,j,i*order+j, j*order+i,dw[i*order+j],dw[j*order+i]);
    }
  //check precomputed derivatives (manual)
  // for(int i=0;i<order;i++)
  // 	printf("prod[%d]=%lf\n",i,prods[i]);
  // for(int i=0;i<order;i++)
  // 	printf("nodes[%d]=%lf\n",i,nodesx[i]);
  // for(int i=0;i<order*order;i++)
  // 	printf("dw[%d]=%lf\n",i, dw[i]);
  // for(int i=0;i<order;i++)
  // printf("%d %le %le\n",i, nodes(i), quadratureWeights(i));
}
double Differential::nodes(unsigned short index)
{
  if(index>=order)
    throw(std::range_error("requested point not a valid mesh index"));
  // printf("requested node %d: %le\n",index, nodesx[index]);
  return nodesx[index];
}
double Differential::quadratureWeights(unsigned short index)
{
  if(index>=order)
    throw(std::range_error("requested point not a valid mesh index"));
  // printf("requested weight %d: %le\n",index, qw[index]);
  return qw[index];
}
double Differential::differentiationWeights(unsigned short polynomial, unsigned short point)
{
  if((polynomial>=order)||(point>=order))
    throw(std::range_error("invalid result point or invalid interpolation polynomiales"));
  return dw[polynomial*order+point];
}
