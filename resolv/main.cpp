#include <cstdlib>
#include <cstdio>

#include <gsl/gsl_vector.h>
#include "eigen.h"
#include "interpolation.h"
#include "differential.h"

constexpr size_t sizes=21;
constexpr size_t sizem0=113;
constexpr size_t sizemu=sizem0/5;

constexpr double s[sizes]={5187000,5310000,5582000,5995000,6254000,6512000,7310000,8317000,9802000,8744000,8059000,7556000,8429000,9839000,7713000,6369000,5063000,3196000,1408000,394000,65000};
constexpr double mu0[sizem0]={0.06764,0.00038,0.00024,0.00019,0.00013,0.00010,0.00010,0.00010,0.00009,0.00008,0.00007,0.00007,0.00008,0.00009,0.00012,0.00015,0.00019,0.00024,0.00029,0.00036,0.00042,0.00047,0.00050,0.00053,0.00054,0.00055,0.00055,0.00055,0.00054,0.00055,0.00057,0.00059,0.00061,0.00064,0.00068,0.00072,0.00075,0.00078,0.00082,0.00090,0.00100,0.00110,0.00120,0.00129,0.00140,0.00155,0.00171,0.00190,0.00211,0.00233,0.00255,0.00280,0.00308,0.00339,0.00373,0.00412,0.00454,0.00498,0.00540,0.00585,0.00639,0.00709,0.00795,0.00886,0.00973,0.01070,0.01182,0.01295,0.01415,0.01543,0.01684,0.01846,0.02025,0.02205,0.02388,0.02610,0.02889,0.03233,0.03649,0.04134,0.04680,0.05307,0.06025,0.06839,0.07766,0.08810,0.09941,0.11156,0.12500,0.14002,0.15698,0.17602,0.19751,0.22205,0.24801,0.27055,0.29434,0.31910,0.34485,0.37165,0.39954,0.42855,0.45874,0.49015,0.52284,0.55684,0.59223,0.62905,0.66736,0.70722,0.74869,0.79185,0.83675};

double mu[sizemu];
double xm[sizemu];
double xs[sizes];

Differential::Type inttype=Differential::Type::ClenshawCurtis;
size_t maxsize=100;
size_t step=1;
size_t start=2;
size_t order=0;
bool analytic=false;

double hint=0;

int main(int argc, char**argv)
{
    for(int i=0;i<sizes;i++)
        xs[i]=5*i;
    for (int i = 0; i < sizemu; ++i) {
        xm[i]=5*i;
        mu[i]=mu0[5*i];
    }

    fprintf(stderr,"usage: analytic(1) order start step maxsize type(g/cc) hint\n");
    switch (argc-1){
    default:
        fprintf(stderr,"Ignoring arguments beyond the fourth (%d supplied)\n",argc);
        [[fallthrough]];
    case 7:
        hint=atof(argv[7]);
        fprintf(stderr, "hint: %e ", hint);
        [[fallthrough]];
    case 6:
        switch(atoi(argv[6])){
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
    case 5:
        maxsize=atoi(argv[5]);
        [[fallthrough]];
    case 4:
        step=atoi(argv[4]);
        if(step==0) {
            step=1;
            fprintf(stderr,"Invlid step 0, setting 1\n");
        }
        [[fallthrough]];
    case 3:
        start=atoi(argv[3]);
        [[fallthrough]];
    case 2:
        order=atoi(argv[2]);
        [[fallthrough]];
    case 1:
        if(atoi(argv[1])) {
            analytic=true;
            fprintf(stderr, "using analytic functions\n");
        } else {
            fprintf(stderr, "using data interpolants\n");
        }
        [[fallthrough]];
    case 0:
        ;
    }

    if(order && !analytic) {
        start=((start-1)/20+1)*20;
        step=((step-1)/20+1)*20;
    }
    
    interpolation::splineInit(sizes, sizemu, xm, mu, xs, s, analytic);

    double* mesh;
    if (order)
        mesh=(double*) malloc((maxsize+1)*sizeof(double));
    else
        mesh=NULL;
    
    for(int i=start; i<=maxsize;i+=step){
        if (order) {
            for (int j = 0; j <= i; ++j)
                mesh[j]=j*100./i;
            eigen::init(order+1, inttype, mesh, i+1);
        } else if (analytic) {
            eigen::init<2>(i, inttype, {0,100});
        } else {
            eigen::init(i, inttype, xs, sizes);
        }
        
        auto a=eigen::computeSpectralRadius(hint);
        printf("%4d %.20e\n", i, a);
        fflush(stdout);
    }
    free(mesh);
    eigen::freemem();
    interpolation::free();
}
