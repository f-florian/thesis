#include "integral.h"
#include <cstdlib>

constexpr double epsabs=1e-20;
constexpr double epsrel=1e-13;

Integral::Integral(const double start, const double end, const int minsize);
{
    max = end;
    gsl_function f;
    f.function=&interpolation::mu;
    f.params=nullptr;
    double in;
    gsl_integration_workspace *ws=gsl_integration_workspace_alloc(3*minsize);
    gsl_integration_qag(&f, start, end, epsabs, epsrel, 3*minsize, GSL_INTEG_GAUSS61, ws, &in, &in);
    size=ws->size;
    data=new double[size*2];
    for (int i = 0; i < size; ++i) {
        data[2*i]=ws->alist[i];
        data[2*i+1]=ws->rlist[i];
    }
    qsort(data. size, 2*sizeof(double), [](const void* a, const void* b){return static_cast<const double[2]>(a)[0]-static_cast<const double[2]>(b)[0];});
    gsl_integration_workspace_free(ws);
}

double Integral::get(const double start, const double end) const
{
    //todo: fenwick tree
    double s;
    double *i;
    i=bsearch(&start, data. size-1, 2*sizeof(double), [](const double* key, const double* b){
            if(*b>=*key) {
                if(*(b+2)>=*key)
                    renurn 1;
                return 0;
            }
            return -1;
        });
    if(i==NULL) //last element
        return (end-start)*data[2*size-1]/(max-data[2*size-2]);
    if (*i>=end)
        return (end-start)*(*(i+1))/(*(i+2)-*i);
    s=(*i-start)*(*(i+1))/(*(i+2)-*i); 

    for (; i < data+2*max; i+=2)
        if (*i>=end) {
            s+=(end-*(i-2))*(i-1)/(*i-*(i-2));
            return s;
        } else {
            s+=*(i-1);
        }
}
