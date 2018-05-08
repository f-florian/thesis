#include <cstdlib>
#include <cstdio>
#include <gsl/gsl_vector.h>
#include "eigen.h"
#include "differential.h"
#include "functions.h"

Differential::Type inttype=Differential::Type::ClenshawCurtis;
size_t maxsize=100;
size_t step=1;
size_t start=2;
size_t order=0;
double length=1;
double c1=0;
double D1=0;

int main(int argc, char**argv)
{
    fprintf(stderr,"usage: order start step maxsize type(g/cc) length D_1 c_1\n");
    switch (argc-1){
    default:
        fprintf(stderr,"Too many arguments (%d supplied); ignoring the last ones\n",argc);
        [[fallthrough]];
    case 8:
        c1=atof(argv[8]);
        fprintf(stderr, "c_1: %e ", c1);
        [[fallthrough]];
    case 7:
        D1=atof(argv[7]);
        fprintf(stderr, "D_1: %e ", D1);
        [[fallthrough]];
    case 6:
        length=atof(argv[6]);
        fprintf(stderr, "length: %e ", length);
        [[fallthrough]];
    case 5:
        switch(atoi(argv[5])){
        case 0:
            inttype=Differential::Type::Gauss;
            fprintf(stderr, "Integration type: Gauss; ");
            break;
        case 1:
            inttype=Differential::Type::ClenshawCurtis;
            fprintf(stderr, "Integration type: Clenshaw-Curtis; ");
            break;
        default:
            fprintf(stderr,"Invalid integration type specified (%d), using default 'Clenshaw-Curtis (1)'\n", atoi(argv[6]));
        }
        [[fallthrough]];
    case 4:
        maxsize=atoi(argv[4]);
        [[fallthrough]];
    case 3:
        step=atoi(argv[3]);
        if(step==0) {
            step=1;
            fprintf(stderr,"Invlid step 0, setting 1\n");
        }
        [[fallthrough]];
    case 2:
        start=atoi(argv[2]);
        if (start>maxsize)
            fprintf(stderr,"minimum size (%ld) > maximum size (%ld); nothing to compute", start, maxsize);
        [[fallthrough]];
    case 1:
        order=atoi(argv[1]);
        [[fallthrough]];
    case 0:
        ;
    }

    //init
    parameters::init(length,c1,D1);
    
    //alloc mesh
    double* mesh;
    if (order)
        mesh=(double*) malloc((maxsize+1)*sizeof(double));
    else
        mesh=NULL;

    //compute
    for(size_t i=start; i<=maxsize;i+=step){
        fprintf(stderr, "%4ld\n", i);
        if (order) {
            for (size_t j = 0; j <= i; ++j)
                mesh[j]=static_cast<double>(j)*100./static_cast<double>(i);
            eigen::init(order+1, inttype, mesh, i+1);
        } else {
            eigen::init<2>(i, inttype, {0, length});
        }
        
        auto a=eigen::computeSpectralRadius(0);
        printf("%4ld %.20e\n", i, a);
        fflush(stdout);
    }
    free(mesh);
    eigen::freemem();
}
