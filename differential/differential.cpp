#include <cmath>
#include <stdexcept>
#include <gsl/gsl_integration.h>
#include "differential.h"

constexpr double pi=3.14159265358979323846264338327950288419716939937510;

Differential::Differential(unsigned short order_p)
{
    //set order
    order=order_p;
    //allocate matrix of derivative weights
    dw=new double[order*order];
    //get quadrature weights from gsl (swap pointers to avoid copying, then delete gsl objects)
    gsl_integration_fixed_workspace *iws=gsl_integration_fixed_alloc(gsl_integration_fixed_legendre, order, 0,1,0,0);
    nodesx=iws->x;
    iws->x=NULL;
    qw=iws->weights;
    iws->weights=NULL;
    gsl_integration_fixed_free(iws);
#ifdef ___DEBUG_PRINT
    for(int i=0;i<order;i++)
        printf("%d %le %le\n",i, nodesx[i], qw[i]);
#endif // ___DEBUG_PRINT
    // store temporary values to accelerate computation:
    double difftmp[order*order];
    double prods[order];                //weights for the barycentric lagrange interpolation
    for(int i=0;i<order;i++){
        prods[i]=1;
        dw[i*order+i]=0;
        // store differences in the lower triagular part; todo: only alloc lower triangular
        for(int j=0;j<i;j++){
            // difftmp_{i,j}=x_i-x_j
            difftmp[i*order+j]=nodesx[i]-nodesx[j];
            // prods_i=\sum_k 1/(x_i-x_k)
            prods[i]*=difftmp[i*order+j];
            prods[j]*=-difftmp[i*order+j];
        }
    }
#ifdef ___DEBUG_PRINT
    for(int i=0;i<order*order;i++)
        printf("dw[%d]=%lf\n",i, dw[i]);
#endif // ___DEBUG_PRINT
    for(int i=0;i<order;i++)
        for(int j=0;j<i;j++) {
            // li'(xj)=TODO write formula (computed values are ok)
            dw[i*order+j]=-prods[j]/prods[i]/difftmp[i*order+j];
            // lj'(xi)=TODO write formula (computed values are ok)
            dw[j*order+i]=prods[i]/prods[j]/difftmp[i*order+j];
#ifdef ___DEBUG_PRINT
            printf("set (%d,%d) (=%d,%d) to %lf (and %lf)\n",i,j,i*order+j, j*order+i,dw[i*order+j],dw[j*order+i]);
#endif // ___DEBUG_PRINT
        }
    for (int i = 0; i < order; ++i) {
        dw[i*order+i]=0;
        for (int j = 0; j<i; ++j) {
            dw[i*order+i]-=dw[j*order+i];
            dw[j*order+j]-=dw[i*order+j];
        }
    }


#ifdef ___DEBUG_PRINT
    for(int i=0;i<order;i++)
        printf("prod[%d]=%lf\n",i,prods[i]);
    for(int i=0;i<order;i++)
        printf("nodes[%d]=%lf\n",i,nodesx[i]);
    for(int i=0;i<order*order;i++)
        printf("dw[%d]=%lf\n",i, dw[i]);
    for(int i=0;i<order;i++)
        printf("qw(%d) %le %le\n",i, nodes(i), quadratureWeights(i));
#endif // ___DEBUG_PRINT
}
double Differential::nodes(unsigned short index, double end, double start)
{
    if(index>=order)
        throw(std::range_error("requested point not a valid mesh index"));
#ifdef ___DEBUG_PRINT
    printf("requested node %d: %le\n",index, nodesx[index]);
#endif // ___DEBUG_PRINT
    // scaling nodes trough affinity
    return start+(end-start)*nodesx[index];
}
double Differential::quadratureWeights(unsigned short index, double end, double start)
{
    if(index>=order)
        throw(std::range_error("requested point not a valid mesh index"));
#ifdef ___DEBUG_PRINT
    printf("requested qweight %d: %le\n",index, qw[index]);
#endif // ___DEBUG_PRINT
    // scaling: \int_a^b f(x) dx=(b-a)\int_0^1 f(a+(b-a)x)dx; nodes are as above
    return (end-start)*qw[index];
}
double Differential::differentiationWeights(unsigned short polynomial, unsigned short point, double end, double start)
{
    if((polynomial>=order)||(point>=order))
        throw(std::range_error("invalid result point or invalid interpolation polynomiales"));
#ifdef ___DEBUG_PRINT
    printf("requested dweight %d,%d: %le\n",polynomial,point, dw[polynomial*order+point]);
#endif // ___DEBUG_PRINT
    // scaling: 
    return dw[polynomial*order+point];
}
double* Differential::StealNodes()
{
    auto tmp=nodesx;
    nodesx=nullptr;
    return tmp;
}
double*  Differential::StealQuadratureWeights()
{
    auto tmp=qw;
    qw=nullptr;
    return tmp;
}
double*  Differential::StealDifferentiationWeights()
{
    auto tmp=dw;
    dw=nullptr;
    return tmp;
}
