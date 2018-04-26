#include <cstdlib>
#include <cstdio>

#include <gsl/gsl_vector.h>
#include "eigen.h"
#include "interpolation.h"
#include "differential.h"

constexpr size_t sizes=21;
constexpr size_t sizem0=113;

constexpr double s[sizes]={5187000,5310000,5582000,5995000,6254000,6512000,7310000,8317000,9802000,8744000,8059000,7556000,8429000,9839000,7713000,6369000,5063000,3196000,1408000,394000,65000};
double mu0[sizem0]={0.06764,0.00038,0.00024,0.00019,0.00013,0.00010,0.00010,0.00010,0.00009,0.00008,0.00007,0.00007,0.00008,0.00009,0.00012,0.00015,0.00019,0.00024,0.00029,0.00036,0.00042,0.00047,0.00050,0.00053,0.00054,0.00055,0.00055,0.00055,0.00054,0.00055,0.00057,0.00059,0.00061,0.00064,0.00068,0.00072,0.00075,0.00078,0.00082,0.00090,0.00100,0.00110,0.00120,0.00129,0.00140,0.00155,0.00171,0.00190,0.00211,0.00233,0.00255,0.00280,0.00308,0.00339,0.00373,0.00412,0.00454,0.00498,0.00540,0.00585,0.00639,0.00709,0.00795,0.00886,0.00973,0.01070,0.01182,0.01295,0.01415,0.01543,0.01684,0.01846,0.02025,0.02205,0.02388,0.02610,0.02889,0.03233,0.03649,0.04134,0.04680,0.05307,0.06025,0.06839,0.07766,0.08810,0.09941,0.11156,0.12500,0.14002,0.15698,0.17602,0.19751,0.22205,0.24801,0.27055,0.29434,0.31910,0.34485,0.37165,0.39954,0.42855,0.45874,0.49015,0.52284,0.55684,0.59223,0.62905,0.66736,0.70722,0.74869,0.79185,0.83675};

constexpr size_t sizemu=sizem0/5;
double mu[sizemu];
double xm[sizemu];

constexpr size_t start=2;
constexpr size_t step=1;

Differential::Type inttype=Differential::Type::ClenshawCurtis;
size_t maxsize=100;
bool analytic=false;

int main(int argc, char**argv)
{
    fprintf(stderr,"usage: analytic(1) maxsize type(g/cc)\n");
    double xs[sizes];
    for(int i=0;i<sizes;i++)
        xs[i]=5*i;
    for (int i = 0; i < sizemu; ++i) {
        xm[i]=5*i;
        mu[i]=mu0[5*i];
    }

    switch (argc-1){
    default:
        fprintf(stderr,"Ignoring arguments beyond the fourth (%d supplied)\n",argc);
        [[fallthrough]];
    case 3:
        switch(atoi(argv[3])){
        case 0:
            inttype=Differential::Type::Gauss;
            break;
        case 1:
            inttype=Differential::Type::ClenshawCurtis;
            break;
        default:
            fprintf(stderr,"Invalid integration type specified (%d), using default 'Clenshaw-Curtis (1)'\n", atoi(argv[3]));
        }
        [[fallthrough]];
    case 2:
        maxsize=atoi(argv[2]);
        if(maxsize==0)
            fprintf(stderr,"You selected 0 as maximum size, no matrix radius will be computed\n");
        [[fallthrough]];
    case 1:
        if(atoi(argv[1]))
            analytic=true;
        [[fallthrough]];
    case 0:
        ;
    }
    
    interpolation::splineInit(sizes, sizemu, xm, mu, xs, s, analytic);

    for(int i=start; i<=maxsize;i+=step){
        // size,start,end
        eigen::init(i, 0, 100, inttype, xs, sizes);
        auto a=eigen::computeSpectralRadius();
        printf("%4d %.20e\n", i, a);
    }
    eigen::freemem();
    interpolation::free();
}
